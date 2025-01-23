
"""As the sphinx markdown for read the docs has a horrific way of displaying and formatting tables, we take the
easier to write csv-tables and convert them into the better formatted one which are terrible to write by hand."""

table_file = "/Users/dennis/Dev/FuFis/RST_Tables/TableCSV.txt"
out_file = table_file.replace('.txt', '_formed.txt')

header = [x.strip() for x in open(table_file).read().split(':header: "')[1].split('\n')[0].replace('"', '').split(',')]
rows = [header] + [[y.strip() for y in x.split(',')] for x in open(table_file).read().split('\n\n')[1].split('\n')]
cols = len(rows[0])
longest_entries = [max([len(x[i]) for x in rows]) for i in range(cols)]

vertical_line = '+'
for c in range(cols):
    vertical_line += '-' * longest_entries[c] + '+'

with open(out_file, 'w') as output:
    output.write(vertical_line + '\n')
    for r, row in enumerate(rows):
        for c in range(cols):
            output.write('|' + row[c] + ' ' * (longest_entries[c] - len(row[c])))
        output.write('|\n')
        if r == 0:
            output.write(vertical_line.replace('-', '='))
        else:
            output.write(vertical_line)
        output.write('\n')


