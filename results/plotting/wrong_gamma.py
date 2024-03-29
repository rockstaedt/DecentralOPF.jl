import matplotlib.pyplot as plt
import pandas as pd

from config import *

for case in ['wrong_weight', 'big_gamma']:

    df_duals = pd.read_csv(
        f'../{case}_duals.csv',
        delimiter=DELIMITER,
        na_filter=False
    )

    iterations = df_duals.iteration.unique()

    fig, ax = plt.subplots()

    lambdas = df_duals.loc[
        (df_duals.dual == 'lambda') & (df_duals.timestep == 1),
        'value'
    ].to_list()

    ax.plot(iterations, lambdas)

    ax.set_ylabel('Lambda for t=1')
    ax.set_xlabel('Iterations')
    ax.grid(True)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR + f'{case}_lambda_t1.png', dpi=DPI)

    if case == 'big_gamma':
        fig, ax = plt.subplots()

        ax.plot(iterations[50:70], lambdas[50:70])

        ax.set_ylabel('Lambda for t=1')
        ax.set_xlabel('Iterations')
        ax.grid(True)

        plt.tight_layout()
        fig.savefig(OUTPUT_DIR + f'{case}_lambda_t1_i_50-70.png', dpi=DPI)

    if case == 'wrong_weight':
        fig, ax = plt.subplots()

        rhos = df_duals.loc[
            (
                (df_duals.dual == 'rho')
                & (df_duals.timestep == 1)
                & (df_duals.line == '2')
            ),
            'value'
        ].to_list()

        ax.plot(iterations, rhos)

        ax.set_ylabel('Rho for t=1 and l=2')
        ax.set_xlabel('Iterations')
        ax.grid(True)

        plt.tight_layout()
        fig.savefig(OUTPUT_DIR + f'{case}_rho_t1.png', dpi=DPI)