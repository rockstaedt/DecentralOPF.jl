import matplotlib.pyplot as plt
import pandas as pd

from config import *

###################################################################
# Data
###################################################################

df_duals = pd.read_csv('../TNS_duals.csv', delimiter=DELIMITER, na_filter=False)

iterations = df_duals.iteration.unique()

###################################################################
# Lambda
###################################################################

fig, ax = plt.subplots()

lambdas = df_duals.loc[
    (df_duals.dual == 'lambda') & (df_duals.timestep == 1),
    'value'
].to_list()

ax.plot(iterations, lambdas)

ax.set_ylabel('Lambda')
ax.set_xlabel('Iterations')
ax.grid(True)

plt.tight_layout()
fig.savefig(OUTPUT_DIR + 'TNS_lambda_t1.png', dpi=DPI)

###################################################################
# Duals Flow
###################################################################

for dual_flow in ['mue', 'rho']:

    fig, ax = plt.subplots()

    duals = df_duals.loc[
        (df_duals.dual == dual_flow) & (df_duals.timestep == 1),
        :
    ]

    lines = duals.line.unique()

    for l in lines:
        ax.plot(
            iterations,
            duals.loc[duals.line == l, 'value'].to_list(),
            label=f'L{l}'
        )

    ax.set_ylabel(dual_flow.capitalize())
    ax.set_xlabel('Iterations')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR + f'TNS_{dual_flow}_t1.png', dpi=DPI)
