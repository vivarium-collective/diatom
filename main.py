import os
import numpy as np
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

# list of exchange reactions with phosphorus
phosphorus_exchanges = [
    'EX_pi_e',
    # 'EX_pppi_e',
]

remove_exchanges = [
    'EX_photon410_e',
    'EX_photon430_e',
    'EX_photon450_e',
    'EX_photon470_e',
    'EX_photon490_e',
    'EX_photon510_e',
    'EX_photon530_e',
    'EX_photon550_e',
    'EX_photon570_e',
    'EX_photon590_e',
    'EX_photon610_e',
    'EX_photon630_e',
    'EX_photon650_e',
    'EX_photon670_e',
    'EX_photon690_e',
]


def get_redfield_exchanges(model):
    """TODO"""
    nitrogen_exchanges = []
    phosphorus_exchanges = []
    carbon_exchanges = []
    # go through all the exchange and extract the ones which include nitrogen, carbon, phosphorus
    for exchange in model.exchanges:
        pass
        # get the coefficients
        # coeff = exchange.
        #
        # if 'N' in coeff:
        #     nitrogen_exchanges.append(exchange)

    return nitrogen_exchanges, phosphorus_exchanges, carbon_exchanges


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
        phosphorus_change=1.0,
        outdir='out'
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
    for exchange in phosphorus_exchanges:
        ex = model.reactions.get_by_id(exchange)
        ex.lower_bound *= phosphorus_change

    # run file in new conditions
    solution_after = model.optimize()

    # get fluxes and remove 0s
    fluxes_before = solution_before.fluxes
    # fluxes_before = fluxes_before[fluxes_before != 0]
    fluxes_after = solution_after.fluxes
    # fluxes_after = fluxes_after[fluxes_after != 0]

    # get the values of just the exchanges
    shown_exchanges = [ex for ex in exchange_ids if ex not in remove_exchanges]
    exchanges_before = [fluxes_before[ex] for ex in shown_exchanges]
    exchanges_after = [fluxes_after[ex] for ex in shown_exchanges]

    # fig = plt.figure(figsize=(10, 8))
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Define the positions for the bars
    bar_width = 0.35  # Width of each bar
    r1 = np.arange(len(exchanges_before))  # Positions for the first set of bars
    r2 = [x + bar_width for x in r1]  # Positions for the second set of bars

    # Plot the bars
    plt.bar(r1, exchanges_before, width=bar_width, label='before', alpha=0.6)
    plt.bar(r2, exchanges_after, width=bar_width, label='after', alpha=0.6)

    # Add labels and legend
    title = f'c{carbon_change}_p{phosphorus_change}_n{nitrogen_change}'
    # plt.yscale('log')
    plt.title(title)
    plt.xticks(fontsize=5)
    plt.legend()
    plt.xticks([r + bar_width / 2 for r in range(len(shown_exchanges))], shown_exchanges, rotation=90)
    plt.ylabel('Flux Value')

    # use grid for vertical lines
    ax.grid(axis='x')

    # make sure outdir exists
    os.makedirs(outdir, exist_ok=True)

    # Save the figure
    plt.savefig(f'{outdir}/fluxes_comparison_{title}.png', dpi=600)

    plt.show()


def scan_redfield(
        c_change_values=None,
        n_change_values=None,
        p_change_values=None,
        outdir='out'
):
    if p_change_values is None:
        p_change_values = [1]
    if n_change_values is None:
        n_change_values = [1]
    if c_change_values is None:
        c_change_values = [1]

    # perform the scan
    for c in c_change_values:
        for n in n_change_values:
            for p in p_change_values:
                compare_constraints_to_baseline(
                    model_file='models/LL.new.xml',
                    carbon_change=c,
                    nitrogen_change=n,
                    phosphorus_change=p,
                    outdir=outdir
                )

if __name__ == '__main__':
    # compare_constraints_to_baseline(
    #     model_file='models/LL.new.xml',
    #     carbon_change=0.0,
    #     nitrogen_change=2.0,
    #     phosphorus_change=1,
    # )
    scan_redfield(
        c_change_values=[i for i in range(0, 20, 5)],
        # n_change_values=[0, 0.4, 0.8, 1.0, 1.4, 1.8, 2.0],
        p_change_values=[i for i in range(-10, 10, 5)],
        outdir='out_test4'
    )

