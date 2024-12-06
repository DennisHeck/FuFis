============
Test
============



.. autofunction:: BasicPlotter.basic_hist


Below Basicplotter.basichist

.. code-block:: python

   import pandas as pd
   testcode

.. code-block:: python

    avg_flipper_length = pd.DataFrame(penguin_df.groupby('species')['flipper_length_mm'].mean())
BasicPlotter.basic_bars(penguin_df, x_col='species', y_col='flipper_length_mm',
                        x_order=['Chinstrap', 'Adelie', 'Gentoo'], title='Group of pinguins on land is called a waddle',
                        output_path=out_dir, y_label='Flipper length [mm]', rotation=None, palette='glasbey_cool')



%.. automodule:: BasicPlotter
5   :members:
%   :undoc-members:
%  :show-inheritance:

Lowest
