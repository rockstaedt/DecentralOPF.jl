from lib2to3.pgen2.pgen import generate_grammar
import matplotlib.pyplot as plt
import pandas as pd

from config import *

###################################################################
# Data
###################################################################

df_generators = pd.read_csv(
    '../TNS_generators.csv',
    delimiter=DELIMITER,
    na_filter=False
)

iterations = df_generators.iteration.unique()
generators = df_generators.generator.unique()

###################################################################
# Generator
###################################################################


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

for ax, generator in zip([ax1, ax2, ax3, ax4], generators):
    generation_values = df_generators.loc[
        (df_generators.generator == generator) & (df_generators.timestep == 1),
        'generation'
    ].to_list()

    ax.plot(iterations, generation_values)

    ax.set_ylabel('Generation')
    ax.set_xlabel('Iterations')
    if generator == 'pv':
        ax.set_title('PV')
    else:
        ax.set_title(generator.capitalize())
    ax.grid(True)

plt.tight_layout()
fig.savefig(OUTPUT_DIR + f'TNS_generation_t1.png', dpi=DPI)