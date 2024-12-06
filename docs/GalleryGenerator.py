

"""Gather the code that creates plots and output to be shown in the documentation. The code blocks also serve
as template for the code in the documentation itself. This is why the imports are above the clode blocks directly."""

# Block that has to be executed for all.
import src.BasicPlotter as BasicPlotter
import pandas as pd
import seaborn as sns
out_dir = 'docs/gallery/'
penguin_df = sns.load_dataset('penguins')

# Barplot with one bar per group.
# ***BasicPlotter.basic_bars
avg_flipper_length = pd.DataFrame(penguin_df.groupby('species')['flipper_length_mm'].mean())
BasicPlotter.basic_bars(penguin_df, x_col='species', y_col='flipper_length_mm', formats=['png'],
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Group of pinguins on land is called a waddle',
                        output_path=out_dir, y_label='Flipper length [mm]', rotation=None, palette='glasbey_cool')
# ---

# ***BasicPlotter.basic_bars2
avg_bill_length = pd.DataFrame(penguin_df.groupby('species')['bill_length_mm'].mean())
BasicPlotter.basic_bars(avg_bill_length, x_col='species', y_col='bill_length_mm',
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Group of pinguins on land is called a waddle',
                        output_path=out_dir, y_label='Bill length [mm]', rotation=None, palette='glasbey_cool')
# ---

# Other stuff



