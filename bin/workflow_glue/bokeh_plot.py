"""Script from plannotate edited to work with bokeh 3."""
from math import pi

from bokeh.models import ColumnDataSource, HoverTool, Range1d, WheelZoomTool
from bokeh.models.annotations import Label
from bokeh.plotting import figure
import numpy as np
import pandas as pd
from plannotate.bokeh_plot import calc_glyphs, calc_level, calc_num_markers
import plannotate.resources as rsc


def get_bokeh(df, linear=False):
    """Get bokeh from plannotate updated to use Bokeh v3 to match ezcharts."""
    # df = df.fillna("")
    x = 0
    y = 0
    baseradius = .18
    tooltips = """
            <font size="4"><b>@Feature</b> — @Type   @pi_permatch_int</font><br>
            @Description
    """  # noqa: E501
    hover = HoverTool(tooltips=tooltips)
    plotsize = .35
    plotdimen = 800
    x_range = Range1d(-plotsize, plotsize, bounds=(-.5, .5), min_interval=.1)
    y_range = Range1d(-plotsize, plotsize, bounds=(-.5, .5), min_interval=.1)
    toolbar = "right"
    p = figure(
        height=plotdimen, width=plotdimen, title="",
        toolbar_location=toolbar, toolbar_sticky=False,
        match_aspect=True,
        sizing_mode='scale_width', tools=['save', 'pan'],
        x_range=x_range, y_range=y_range)

    # x_range=(-plotsize, plotsize), y_range=(-plotsize, plotsize))
    p.toolbar.logo = None
    p.add_tools(WheelZoomTool(zoom_on_axis=False))
    p.toolbar.active_scroll = p.select_one(WheelZoomTool)

    # backbone line
    p.circle(
        x=x, y=y, radius=baseradius, line_color="#000000",
        fill_color=None, line_width=2.5)

    df = calc_level(df)

    if linear:
        line_length = baseradius / 5
        p.line(
            [0, 0], [baseradius - line_length, baseradius + line_length],
            line_width=4, level="overlay", line_color="black")

    df['pi_permatch_int'] = df['pi_permatch'].astype('int')

    df['pi_permatch_int'] = df['pi_permatch_int'].astype(str) + "%"
    # removes percent from infernal hits
    df.loc[df['db'] == "Rfam", 'pi_permatch_int'] = ""

    df['rstart'] = ((df["qstart"]/df["qlen"])*2*pi)
    df['rend'] = ((df["qend"]/df["qlen"])*2*pi)
    df['rstart'] = np.where(
        df['rstart'] < 0, df['rstart'] + (2*pi), df['rstart'])
    df['rend'] = np.where(df['rend'] < 0, df['rend'] + (2*pi), df['rend'])
    df['rend'] = np.where(
        df['rend'] < df['rstart'], df['rend'] + (2*pi), df['rend'])

    df['Type'] = df['Type'].str.replace('rep_origin', 'origin of replication')

    # DDE0BD
    # C97064
    # C9E4CA
    fullcolordf = pd.read_csv(
        rsc.get_resource("data", "colors.csv"), index_col=0)
    fragcolordf = fullcolordf.copy()
    fragcolordf[['fill_color', 'line_color']] = fragcolordf[
        ['line_color', 'fill_color']]
    fragcolordf["fill_color"] = "#ffffff"

    full = df[df["fragment"] == False] # noqa

    full = full.merge(
        fullcolordf,
        how="left", on=["Type"])
    full['legend'] = full['Type']
    full = full.fillna(
        {"color": "grey",
         "fill_color": "#808080",
         "line_color": "#000000"})
    frag = df[df["fragment"]]
    frag = frag.merge(fragcolordf, how="left", on=["Type"])
    frag = frag.fillna(
        {"color": "grey",
         "fill_color": "#ffffff",
         "line_color": "#808080"})

    df = full.append(frag).reset_index(drop=True)

    # add orientation column
    orient = pd.read_csv(
        rsc.get_resource("data", "feature_orientation.csv"),
        header=None, names=["Type", "has_orientation"])
    orient['Type'] = orient['Type']
    orient['has_orientation'] = orient['has_orientation'].map({"T": True})
    df = df.merge(orient, on="Type", how="left")
    df['Type'] = df['Type'].str.replace("_", " ")
    df['has_orientation'] = df['has_orientation'].fillna(value=False)
    df[[
        'x', 'y', "Lx1", "Ly1",
        "annoLineColor", "lineX",
        "lineY", "theta", "text_align"]] = df.apply(calc_glyphs, axis=1)
    df['legend'] = df['Type']
    # allowedtypes = ['CDS',"promoter","origin of replication","swissprot"]
    allowedtypes = fullcolordf['Type']
    mask = ~df['legend'].isin(allowedtypes)
    df.loc[mask, 'legend'] = 'misc feature'

    # plot annotations
    source = ColumnDataSource(df)
    hover_labels = p.patches(
        'x', 'y', fill_color='fill_color', line_color='line_color',
        name="features", line_width=2.5,
        source=source, legend_group="legend")
    p.multi_line(
        xs="lineX", ys="lineY", line_color="annoLineColor",
        line_width=3, level="overlay",
        line_cap='round', alpha=.5, source=source)

    # `text_align` cannot read from `source` -- have to do this workaround
    right = ColumnDataSource(df[df['text_align'] == 'right'])
    left = ColumnDataSource(df[df['text_align'] == 'left'])
    bcenter = ColumnDataSource(df[df['text_align'] == 'b_center'])
    tcenter = ColumnDataSource(df[df['text_align'] == 't_center'])

    text_level = 'overlay'
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=3, y_offset=8,
        text_align="left",
        text='Feature', level=text_level, source=right)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=-5, y_offset=8,
        text_align="right", text='Feature', level=text_level, source=left)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=0, y_offset=15,
        text_align="center", text='Feature',
        level=text_level, source=bcenter)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=0, y_offset=0,
        text_align="center",
        text='Feature', level=text_level, source=tcenter)

    # calculate chunk size(s) for drawing lines
    plaslen = df.iloc[0]['qlen']
    ticks = calc_num_markers(plaslen)
    ticks_cds = ColumnDataSource(ticks)
    p.multi_line(
        xs="lineX", ys="lineY", line_color="black", line_width=2,
        level="underlay", line_cap='round',
        alpha=.5, source=ticks_cds)

    right = ColumnDataSource(ticks[ticks['text_align'] == 'right'])
    left = ColumnDataSource(ticks[ticks['text_align'] == 'left'])
    bcenter = ColumnDataSource(ticks[ticks['text_align'] == 'b_center'])
    tcenter = ColumnDataSource(ticks[ticks['text_align'] == 't_center'])
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=3, y_offset=6,
        text_align="left", text='bp', alpha=.5,
        text_font_size='size', level=text_level, source=right)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=-5,
        y_offset=6, text_align="right",
        text='bp', alpha=.5, text_font_size='size',
        level=text_level, source=left)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=0, y_offset=15,
        text_align="center", text='bp', alpha=.5,
        text_font_size='size', level=text_level, source=bcenter)
    p.text(
        x="Lx1", y="Ly1", name="2", x_offset=0, y_offset=-3,
        text_align="center", text='bp', alpha=.5,
        text_font_size='size', level=text_level, source=tcenter)
    p.add_tools(hover)
    p.add_layout(Label(
        x=0, y=0, name="2", x_offset=0, y_offset=-8,
        text_align="center",
        text=f"{plaslen} bp", text_color="#7b7b7b",
        text_font_size='16px', level=text_level))
    p.hover.renderers = [hover_labels]
    p.axis.axis_label = None
    p.axis.visible = False
    p.grid.grid_line_color = "#EFEFEF"
    p.outline_line_color = "#DDDDDD"
    p.legend.location = 'bottom_left'
    p.legend.border_line_color = "#EFEFEF"
    p.legend.visible = True
    return p
