import cobra
from cobra import Metabolite, Reaction
import matplotlib.pyplot as plt
import pandas as pd

# list of exchange reactions with nitrogen
nitrogen_exchanges = []

# list of exchange reactions with carbon
carbon_exchanges = [
    'EX_co2_e'
]

# list of exchange reactions with phosphorous
phosphorous_exchanges = [
    'EX_pi_e',
    # 'EX_pppi_e',
]


def get_redfield_exchanges(model):
    """TODO"""
    nitrogen_exchanges = []
    phosphorous_exchanges = []
    carbon_exchanges = []
    # go through all the exchange and extract the ones which include nitrogen, carbon, phosphorous
    for exchange in model.exchanges:
        pass
        # get the coefficients
        # coeff = exchange.
        #
        # if 'N' in coeff:
        #     nitrogen_exchanges.append(exchange)

    return nitrogen_exchanges, phosphorous_exchanges, carbon_exchanges


def update_model_exchanges(model):
    # update model exchanges
    metabolite_to_add_to_e = [
        {'id': 'dgsn', 'compartment': 'c'},
        {'id': '4abz', 'compartment': 'm'},
        {'id': 'ump', 'compartment': 'c'},
        {'id': '4abut', 'compartment': 'c'},
        {'id': '4ahmmp', 'compartment': 'c'},
        {'id': 'adn', 'compartment': 'c'},
        {'id': '4aammp', 'compartment': 'c'},
        {'id': 'arg__L', 'compartment': 'c'},
        {'id': 'chtbs', 'compartment': 'c'},
        {'id': 'Lcyst', 'compartment': 'c'},
        {'id': 'cys__L', 'compartment': 'c'},
        {'id': 'cytd', 'compartment': 'c'},
        {'id': 'g6p', 'compartment': 'c'},
        {'id': 'gthrd', 'compartment': 'c'},
        {'id': 'gsn', 'compartment': 'c'},
        {'id': 'hom__L', 'compartment': 'c'},
        {'id': 'lys__L', 'compartment': 'c'},
        {'id': 'met__L', 'compartment': 'c'},
        {'id': 'phe__L', 'compartment': 'c'},
        {'id': 'ptrc', 'compartment': 'c'},
        {'id': 'sarcs', 'compartment': 'c'},
        {'id': 'glyc3p', 'compartment': 'c'},
        {'id': 'trp__L', 'compartment': 'c'},
        {'id': 'tyr__L', 'compartment': 'c'},
        {'id': 'uri', 'compartment': 'c'}
    ]

    for met in metabolite_to_add_to_e:
        id = met['id']
        compart = met['compartment']

        internal_id = f'{id}_{compart}'
        external_id = f'{id}_e'

        # Create the external metabolite
        external_metabolite = Metabolite(external_id, name=id, compartment='e')

        # Retrieve the internal metabolite
        internal_metabolite = model.metabolites.get_by_id(internal_id)

        # Create the transport reaction
        transport_reaction = Reaction(f'{id}_e_transport')
        transport_reaction.add_metabolites({
            internal_metabolite: -1,
            external_metabolite: 1
        })
        transport_reaction.lower_bound = -1000
        transport_reaction.upper_bound = 1000

        # Add the transport reaction to the model
        model.add_reactions([transport_reaction])

        # Add the exchange reaction for the external metabolite
        exchange_reaction = model.add_boundary(external_metabolite, type="exchange")
        exchange_reaction.lower_bound = 0
        exchange_reaction.upper_bound = 1000

        # Check if the exchange reaction was added and its bounds
        assert exchange_reaction.lower_bound == 0, f"{exchange_reaction.id} has incorrect lower bound."
        assert exchange_reaction.upper_bound == 1000, f"{exchange_reaction.id} has incorrect upper bound."

    return model


def compare_constraints_to_baseline(
        model_file='models/LL.new.xml',
        carbon_change=1.0,
        nitrogen_change=1.0,
        phosphorous_change=1.0
):

    # load model file
    model = cobra.io.read_sbml_model(model_file)

    # # update exchanges
    # model = update_model_exchanges(model)
    exchange_ids = [ex.id for ex in model.exchanges]

    # run file in default conditions
    solution_before = model.optimize()

    # update constraints
    # TODO -- make this more general to any desired exchanges
    for exchange in carbon_exchanges:
        ex = model.reactions.get_by_id(exchange)
        ex.lower_bound *= carbon_change
    for exchange in nitrogen_exchanges:
        ex = model.reactions.get_by_id(exchange)
        ex.lower_bound *= nitrogen_change
    for exchange in phosphorous_exchanges:
        ex = model.reactions.get_by_id(exchange)
        ex.lower_bound *= phosphorous_change

    # run file in new conditions
    solution_after = model.optimize()

    # get fluxes and remove 0s
    fluxes_before = solution_before.fluxes
    # fluxes_before = fluxes_before[fluxes_before != 0]
    fluxes_after = solution_after.fluxes
    # fluxes_after = fluxes_after[fluxes_after != 0]

    # Define the width of the bars
    bar_width = 0.4

    # make a pandas bar graph from the fluxes and save to file
    # plt.style.use('dark_background')
    fig = plt.figure()
    # plot_before = [id for fluxes_before.index if id in exchange_ids]
    values_before = [fluxes_before[ex] for ex in exchange_ids]
    values_after = [fluxes_after[ex] for ex in exchange_ids]
    plt.bar(exchange_ids, values_before, label='before', alpha=0.6)
    plt.bar(exchange_ids, values_after, label='after', alpha=0.6)
    plt.legend()
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.ylabel('Flux Value')
    plt.savefig(f'fluxes_comparison_c{carbon_change}_p{phosphorous_change}_n{nitrogen_change}.png', dpi=300)


if __name__ == '__main__':
    compare_constraints_to_baseline(
        model_file='models/LL.new.xml',
        carbon_change=-0.1,
        nitrogen_change=1,
        phosphorous_change=1,
    )
