import pandas as pd
import plotly.graph_objects as go
import random

with open(r'./ferro_search_results/results.info', encoding='utf-8') as f:
    info = []
    info_block = []
    for line in f:
        if line.startswith('-') and info_block == []:
            continue
        elif line.startswith('-') and info_block != []:
            info.append(info_block)
            info_block = []
        else:
            info_block.append(line)

energy_info = []
max_p_info = []
spg_info = []
id_info = []

for i in info:
    value = []
    p = ''
    id = ''
    if i[1].startswith('Polar_value'):
        energy_info.append(float(i[0].split(' ')[5]))
        id = i[0].rstrip('\n').split(' ')[1]
        p = i[1].split(' ')
        value.append(p[1])
        value.append(p[2])
        value.append(p[3])
    else:
        continue
    index = value.index(max(value))
    lattice_parameters = []
    lattice_info = i[6].rstrip('\n').split(' ')
    lattice_parameters.append(lattice_info[1])
    lattice_parameters.append(lattice_info[2])
    lattice_parameters.append(lattice_info[3])
    polar_lattice = lattice_parameters[index]
    polar = round(abs(float(max(value)) / float(polar_lattice)), 2)
    max_p_info.append(polar)
    spg = i[5].rstrip('\n').split(' ')[2]
    spg_info.append(spg)
    id_info.append(id)
data = {'energy': energy_info, 'polar': max_p_info, 'spg': spg_info, 'id': id_info}
data = pd.DataFrame(data)


unique_spgs = data['spg'].unique()
random_colors = ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(len(unique_spgs))]
spg_color_map = dict(zip(unique_spgs, random_colors))
data['color'] = data['spg'].map(spg_color_map)
fig = go.Figure()
for spg, color in spg_color_map.items():
    spg_data = data[data['spg'] == spg]
    fig.add_trace(go.Scatter(
        x=spg_data['energy'],
        y=spg_data['polar'],
        mode='markers',
        name=f'SPG: {spg}',
        marker=dict(size=12, color=color, line=dict(width=2)),
        hovertemplate='Energy: %{x}<br>Polar: %{y}<br>SPG: %{text}<br>ID: %{customdata}<extra></extra>',
        text=spg_data['spg'],
        customdata=spg_data['id']
    ))
fig.update_layout(
    template='plotly_white',
    font=dict(family='Times New Roman', size=18),
    title_font=dict(size=24, family='Times New Roman', color='black'),
    title='Energy vs Polarization Scatter Plot',
    xaxis=dict(
        title='Energy per Atom (eV)',
        title_font=dict(size=20),
        tickfont=dict(size=16),
    ),
    yaxis=dict(
        title='Abs. Polarization (quantum)',
        title_font=dict(size=20),
        tickfont=dict(size=16)
    ),
    hovermode='closest'
)
fig.add_hline(y=0.0, line_dash='dot', line_color='gray', line_width=2)
fig.add_hline(y=1.0, line_dash='dot', line_color='gray', line_width=2)
fig.write_html(r'./ferro_search_results/scatter_plot.html')
fig.write_image(r'./ferro_search_results/scatter_plot.png')
