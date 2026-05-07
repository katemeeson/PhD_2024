def integrate_omics(model_file, biomass_id, objective_function, growth_threshold, omics, media=None):
    import cobra
    from cobra.io import read_sbml_model
    from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
    from cobra import Model as CobraModel
    
    '''
    INPUT: 
    model_file: file string for COBRA GEM; must be in .xml SBML format; e.g. 'GEM.xml'
    media: dictionary; exchange reaction IDs as keys and g/L as values; default is empty and to use inbuilt medium conditions 
    biomass_id: ID for the biomass equation in model_file; this is unique to each model so must be User specified
    objective_function: chosen objective function; dictionary format; e.g. {model.reactions.get_by_id('biomass_id'):1}
    growth_threshold: a numeric string of growth rate; derived from experimental measurements; units g/gDW/hour which is equivalent to cell doublings per hour
    omics: Pandas DataFrame with sample ID as column and gene/protein IDs as rows; gene or protein IDs must match model IDs

    METHOD: Uses essential reaction and gene lists and experimental growth thresholds to integrate omics data as enzyme expression values (reaction boundaries) according to GPR logic 

    OUTPUT: 
    model: CobraModel; omics-constrained model
    final_fluxes: Pandas DataFrame with fluxes column and Reaction IDs as rows
    constraints_dictionary: 
    '''

    # Step 1: Load COBRA model
    model = read_sbml_model(model_file)
    if not isinstance(model, CobraModel):
        raise TypeError('Model must be a cobra.Model instance') # Check model is in the correct format (CobraModel)

    # Step 2: Input the media conditions and objective function and solve initial FBA
    # Did the User specify a medium and was it in the correct format?
    if media == None:
        media = model.medium
        print('🧪 User did not specify medium. So, the in-built model medium conditions are being used.')
        print('🧪 These are:')
        for k,v in model.medium.items():
            print(k, '|', model.reactions.get_by_id(k).name, '|', v,'g/L')
    elif media != None: 
        if isinstance(media, dict): # Check media is in the correct format (dict)
            print('🧪 User has specified medium:')
            for k,v in media.items():
                print(k, '|', model.reactions.get_by_id(k).name, '|', v,'g/L')
        if not isinstance(media, dict):
            raise TypeError('Media should be a dict instance')
    # If specified, then apply the media conditions:
    if media is not None: 
        model.medium = media
        # Check objective function is in the correct format
        if isinstance(objective_function, dict):
            print('The User specified objective function is:',objective_function)
        if not isinstance(objective_function, dict):
            raise TypeError('Objective function should be dict instance')
    # Apply objective function and solve FBA with media constraints only
    model.objective = objective_function
    media_only_simulation = model.optimize()
    # Check that model still grows
    print('Growth rate with media only is:', media_only_simulation.fluxes[biomass_id],'g/gDW/hour')
    if media_only_simulation.fluxes[biomass_id] <= 0:
        print('Model does not grow under specified media conditions, reconsider media?')

    # Calculate essential reactions on this media (list to be used later)
    # This step can take a while depending on model size but is worth it to avoid errors later on
    print('')
    print('Running essential reactions simulations...')
    result = single_reaction_deletion(model)
    media_essential_reactions = [
        list(x)[0]
        for x in result.loc[result['growth'] <= 1e-6, 'ids']
    ]
    print('There are',len(media_essential_reactions),'essential reactions with this medium')

    print('Running essential genes simulations...')
    result_genes = single_gene_deletion(model)
    media_essential_genes = [
        next(iter(x))
        for x in result_genes.loc[result_genes['growth'] <= 1e-6, 'ids']
    ]
    print('There are',len(media_essential_genes),'essential genes with this medium')
    print('')
    
    # Step 3: Provide initial separation of model according to reaction rules
    # Initial separation
    print('')
    ANDs = []
    ORs = []
    ANDORs = []
    one_gene = []
    no_gene = []
    for r in model.reactions:
        if 'and' in r.gene_reaction_rule and 'or' not in r.gene_reaction_rule:
            ANDs.append(r.id)
        if 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule:
            ANDORs.append(r.id)
        if 'or' in r.gene_reaction_rule and 'and' not in r.gene_reaction_rule:
            ORs.append(r.id)
        if len(r.gene_reaction_rule) == 0:
            no_gene.append(r.id)
        if len(r.gene_reaction_rule) != 0:
            if 'or' in r.gene_reaction_rule:
                continue
            elif 'and' in r.gene_reaction_rule:
                continue
            else:
                one_gene.append(r.id)
    print('AND rules: ', len(ANDs))
    print('ANDOR rules: ', len(ANDORs))
    print('OR rules: ', len(ORs))
    print('ONE GENE rules: ', len(one_gene))
    print('NO GENE rules: ', len(no_gene))
    print('Proportion of model not annotated: ', len(no_gene)/len(model.reactions)*100, '%')
    print('Proportion of model which IS annotated: ', 100-(len(no_gene)/len(model.reactions)*100), '%')
    print('Total Reactions = ', len(model.reactions))
    if len(model.reactions) == (len(ORs)+len(ANDs)+len(ANDORs)+len(one_gene)+len(no_gene)):
        print('Model has been completely, successfuly separated')
    else:
        print('Some reactions have not been grouped')
    # Separation according to reversibility
    #separate model further, according to reversibility of rule
    one_gene_forward = []
    one_gene_reversible = []
    for r in one_gene:
        if model.reactions.get_by_id(r).reversibility == True:
            one_gene_reversible.append(r)
        else:
            one_gene_forward.append(r)
    print('Number of reversible, ONE GENE reactions =', len(one_gene_reversible))
    print('Number of forward, ONE GENE reactions =', len(one_gene_forward))
    or_forward = []
    or_reversible = []
    for r in ORs:
        if model.reactions.get_by_id(r).reversibility == True:
            or_reversible.append(r)
        else:
            or_forward.append(r)
    print('Number of reversible, OR reactions =', len(or_reversible))
    print('Number of forward, OR reactions =', len(or_forward))
    and_forward = []
    and_reversible = []
    for r in ANDs:
        if model.reactions.get_by_id(r).reversibility == True:
            and_reversible.append(r)
        else:
            and_forward.append(r)
    print('Number of reversible, AND reactions =', len(and_reversible))
    print('Number of forward, AND reactions =', len(and_forward))
    print('')
    
    # Integrate omics constraints
    # Build enzyme expression dictionary
    expression_dictionary = {
    str(k): float(v)
    for k, v in transcriptomics.iloc[:, 0].to_dict().items() # keys need to be strings so that they match the model gene_reaction_rule type
    }
    
    # Constrain according to reaction rules
    constraints_dictionary = {}
    print('')
    print('ONE GENE forward rules:')
    for r in one_gene_forward:
        if r not in media.keys():
            if r not in media_essential_reactions:
                if model.reactions.get_by_id(r).gene_reaction_rule not in media_essential_genes:
                    if model.reactions.get_by_id(r).gene_reaction_rule in expression_dictionary.keys(): 
                        if expression_dictionary[model.reactions.get_by_id(r).gene_reaction_rule] != (0,0):
                            model.reactions.get_by_id(r).bounds = (0,(expression_dictionary[model.reactions.get_by_id(r).gene_reaction_rule]))
                            solution = model.optimize()
                            if solution.fluxes[biomass_id] <= growth_threshold:
                                print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                model.reactions.get_by_id(r).bounds = (0,1000)
                                solution = model.optimize()
                                print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                            if solution.fluxes[biomass_id] != 0:
                                if solution.fluxes[biomass_id] > growth_threshold:
                                    print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
    print('')
    print('ONE GENE reversible rules:')
    for r in one_gene_reversible:
        if r not in media.keys():
            if r not in media_essential_reactions: 
                if model.reactions.get_by_id(r).gene_reaction_rule not in media_essential_genes:
                    if model.reactions.get_by_id(r).gene_reaction_rule in expression_dictionary.keys(): 
                        if expression_dictionary[model.reactions.get_by_id(r).gene_reaction_rule] != (0,0): 
                            model.reactions.get_by_id(r).bounds = (-1*(expression_dictionary[model.reactions.get_by_id(r).gene_reaction_rule]),(expression_dictionary[model.reactions.get_by_id(r).gene_reaction_rule]))
                            solution = model.optimize()
                            if solution.fluxes[biomass_id] <= growth_threshold:
                                print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                model.reactions.get_by_id(r).bounds = (-1000,1000)
                                solution = model.optimize()
                                print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                            if solution.fluxes[biomass_id] != 0:
                                if solution.fluxes[biomass_id] > growth_threshold:
                                    print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
    print('')
    print('OR forward rules:')
    for r in or_forward:
        if r not in media.keys():
            if r not in media_essential_reactions:
                rule_list = []
                rule = model.reactions.get_by_id(r).gene_reaction_rule
                rule_list.append(rule.split())
                rule_genes = []
                for n in rule_list:
                    for num in n:
                        if num not in ['and', 'or', '(', ')']:
                            rule_genes.append(num.strip('()'))
                #calculate sum of expressions
                genes_in_dataset = []
                for gene in rule_genes:
                    if gene in expression_dictionary.keys():
                        genes_in_dataset.append(gene)
                exp_list = []        
                for g in genes_in_dataset:
                    exp_list.append(expression_dictionary[g])
                sum_of_expressions = sum(exp_list)
                #set bounds according to sum of expressions
                model.reactions.get_by_id(r).bounds = (0,sum_of_expressions)
                solution = model.optimize()
                if solution.fluxes[biomass_id] <= growth_threshold:
                    print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    model.reactions.get_by_id(r).bounds = (0,1000)
                    solution = model.optimize()
                    print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                if solution.fluxes[biomass_id] != 0:
                    if solution.fluxes[biomass_id] > growth_threshold:
                        print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                        constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
    print('')
    print('OR reversible rules:')
    for r in or_reversible:
        if r not in media.keys():
            if r not in media_essential_reactions:
                rule_list = []
                rule = model.reactions.get_by_id(r).gene_reaction_rule
                rule_list.append(rule.split())
                rule_genes = []
                for n in rule_list:
                    for num in n:
                        if num not in ['and', 'or', '(', ')']:
                            rule_genes.append(num.strip('()'))
                #calculate sum of expressions
                genes_in_dataset = []
                for gene in rule_genes:
                    if gene in expression_dictionary.keys():
                        genes_in_dataset.append(gene)
                exp_list = []        
                for g in genes_in_dataset:
                    exp_list.append(expression_dictionary[g])
                sum_of_expressions = sum(exp_list)
                #set bounds according to sum of expressions
                model.reactions.get_by_id(r).bounds = (-1*sum_of_expressions,sum_of_expressions)
                solution = model.optimize()
                if solution.fluxes[biomass_id] <= growth_threshold:
                    print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    model.reactions.get_by_id(r).bounds = (-1000,1000)
                    solution = model.optimize()
                    print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                if solution.fluxes[biomass_id] != 0:
                    if solution.fluxes[biomass_id] > growth_threshold:
                        print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                        constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
    print('')
    print('AND forward rules:')
    for r in and_forward:
        if r not in media.keys():
            if r not in media_essential_reactions:
                rule_list = []
                rule = model.reactions.get_by_id(r).gene_reaction_rule
                rule_list.append(rule.split())
                rule_genes = []
                for n in rule_list:
                    for num in n:
                        if num not in ['and', 'or', '(', ')']:
                            rule_genes.append(num.strip('()'))
                #calculate sum of expressions
                genes_in_dataset = []
                for gene in rule_genes:
                    if gene in expression_dictionary.keys():
                        genes_in_dataset.append(gene)
                exp_list = []        
                for g in genes_in_dataset:
                    exp_list.append(expression_dictionary[g])
                if len(exp_list) > 1:
                    exp_list.sort()
                    min_expression = exp_list[0]
                if len(exp_list) == 1:
                    min_expression = exp_list[0]
                if len(exp_list) == 0:
                    min_expression = 0
                #set bounds according to sum of expressions
                model.reactions.get_by_id(r).bounds = (0,min_expression)
                solution = model.optimize()
                if solution.fluxes[biomass_id] <= growth_threshold:
                    print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    model.reactions.get_by_id(r).bounds = (0,1000)
                    solution = model.optimize()
                    print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                if solution.fluxes[biomass_id] != 0:
                    if solution.fluxes[biomass_id] > growth_threshold:
                        print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                        constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
    print('')
    print('AND reversible rules:')
    for r in and_reversible:
        if r not in media.keys():
            if r not in media_essential_reactions:
                rule_list = []
                rule = model.reactions.get_by_id(r).gene_reaction_rule
                rule_list.append(rule.split())
                rule_genes = []
                for n in rule_list:
                    for num in n:
                        if num not in ['and', 'or', '(', ')']:
                            rule_genes.append(num.strip('()'))
                #calculate sum of expressions
                genes_in_dataset = []
                for gene in rule_genes:
                    if gene in expression_dictionary.keys():
                        genes_in_dataset.append(gene)
                exp_list = []        
                for g in genes_in_dataset:
                    exp_list.append(expression_dictionary[g])
                if len(exp_list) > 1:
                    exp_list.sort()
                    min_expression = exp_list[0]
                if len(exp_list) == 1:
                    min_expression = exp_list[0]
                if len(exp_list) == 0:
                    min_expression = 0
                #set bounds according to sum of expressions
                model.reactions.get_by_id(r).bounds = ((-1*min_expression),min_expression)
                solution = model.optimize()
                if solution.fluxes[biomass_id] <= growth_threshold:
                    print(r, ':', 'constrained bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    model.reactions.get_by_id(r).bounds = (-1000,1000)
                    solution = model.optimize()
                    print(r, ':', 're-opened bounds', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                    constraints_dictionary[r] = model.reactions.get_by_id(r).bounds
                if solution.fluxes[biomass_id] != 0:
                    if solution.fluxes[biomass_id] > growth_threshold:
                        print(r, ':', model.reactions.get_by_id(r).bounds, solution.fluxes[biomass_id])
                        constraints_dictionary[r] = model.reactions.get_by_id(r).bounds

    # Apply the final constraints dictionary to the media-constrained model
    for k,v in constraints_dictionary.items():
        model.reactions.get_by_id(k).bounds = constraints_dictionary[k]
    final_solution = model.optimize()
    print('')
    doubling_time = 1/final_solution.fluxes[biomass_id]
    print('Growth after integrating omics:', doubling_time, 'hours')
    print('equivalent to:', final_solution.fluxes[biomass_id], 'g/gDW/hour')
    final_fluxes = final_solution.fluxes.to_frame()
    
    return(model, final_fluxes, constraints_dictionary)