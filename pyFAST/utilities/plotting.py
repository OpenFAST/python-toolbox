"""Provides the error plotting and html generation functionality."""


import os
from typing import List, Tuple

import numpy as np
from bokeh.embed import components
from bokeh.layouts import gridplot
from bokeh.plotting import ColumnDataSource, figure
from bokeh.models.tools import HoverTool

INDENT = "  "


# Plotting-specific helper functions
def _replace_id_div_string(div: str, attribute: str) -> str:
    """
    Replaces the ID in the stringified div container.

    Parameters
    ----------
    div : str
        Div object as a string.
    attribute : str
        Name of the attribute to be inserted as the ID.

    Returns
    -------
    str
        Updated div string.
    """

    id_start = div.find("id=") + 4
    id_end = div[id_start:].find('"') + id_start
    div = attribute.join((div[:id_start], div[id_end:]))
    # div = INDENT.join((div[:1], div[1:]))
    div = div.replace("\n", f"\n{INDENT}")
    return div


def _replace_id_script_string(script: str, attribute: str) -> str:
    """
    Replaces the ID in the stringified script object with the attribute name.

    Parameters
    ----------
    script : str
        Script object as a string.
    attribute : str
        Name of the attribute the script refers to.

    Returns
    -------
    str
        Script object as a string.
    """
    id_start = script.find("var render_items")
    id_start += script[id_start:].find("roots")
    id_start += script[id_start:].find('":"') + 3
    id_end = script[id_start:].find('"') + id_start
    script = attribute.join((script[:id_start], script[id_end:]))
    script.replace("\n", f"\n{INDENT}")
    return script


def create_header(title: str) -> str:
    """
    Create the html header for the results summary.

    Parameters
    ----------
    title : str
        Title of the document.

    Returns
    -------
    str
        HTML head as a string.
    """

    html = "\n".join(
        (
            "<!DOCTYPE html>",
            "<html>",
            "<head>",
            f"{INDENT}<title>{title}</title>",
            f'{INDENT}<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">',
            f'{INDENT}<link href="https://cdn.pydata.org/bokeh/release/bokeh-widgets-1.4.0.min.css" rel="stylesheet" type="text/css">',
            f'{INDENT}<link href="https://cdn.pydata.org/bokeh/release/bokeh-1.3.4.min.css" rel="stylesheet" type="text/css">',
            f'{INDENT}<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>',
            f'{INDENT}<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>',
            f'{INDENT}<script src="https://cdn.pydata.org/bokeh/release/bokeh-1.4.0.min.js"></script>',
            f'{INDENT}<script src="https://cdn.pydata.org/bokeh/release/bokeh-widgets-1.4.0.min.js"></script>',
            f'{INDENT}<script type="text/javascript"> Bokeh.set_log_level("info"); </script>',
            "",
            f'{INDENT}<style media="screen" type="text/css">',
            f"{INDENT * 2}.cell-warning {{",
            f"{INDENT * 3}background-color: #FF6666;",
            f"{INDENT * 2}}}",
            f"{INDENT * 2}.cell-highlight {{",
            f"{INDENT * 3}background-color: #E5E589;",
            f"{INDENT * 2}}}",
            f"{INDENT}</style>",
            "</head>",
        )
    )
    return html


def create_tail() -> str:
    """Creates the end of the html file."""
    return "</html>\n"


def create_table_head(columns: List[str]) -> str:
    """
    Creates the head for an HTML table.

    Parameters
    ----------
    columns : List[str]
        List of column names.

    Returns
    -------
    str
        Stringified table head that can be inserted into an HTML document.
    """

    header = "\n".join((f"{INDENT * 4}<th>{col}</th>" for col in columns))

    head = "\n".join(
        (
            f'{INDENT * 2}<table class="table table-bordered table-hover table-sm" style="margin: auto; width: 100%; font-size:80%">',
            f"{INDENT * 3}<tr>",
            header,
            f"{INDENT * 3}</tr>",
        )
    )
    return head


# Plotting functions


def plot_single_attribute_error(
    x: np.ndarray, y1: np.ndarray, y2: np.ndarray, xlabel: str, title1: str, title2: str
) -> tuple:
    """
    Plots x vs y1 and x vs y2 as a grid with shared x-axis properties.

    Parameters
    ----------
    x : np.ndarray
        Time, in seconds.
    y1 : np.ndarray
        Baseline solution.
    y2 : np.ndarray
        Test solution.
    xlabel : str
        Label for `x`.
    title1 : str

    title2 : str
        [description]

    Returns
    -------
    script : str
        Text of the script to be embedded in an html file.
    div : str
        Text of the div element to be embedded in an html file.
    """

    # Create the data source
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.absolute(y1 - y2) / y1
    data = ColumnDataSource(data=dict(time=x, baseline=y1, local=y2, error=norm))

    # Plot the unadjusted comparison of values
    p1 = figure(title=title1)
    p1.line(
        "time",
        "baseline",
        color="#1b9e77",
        line_width=3,
        legend_label="Baseline",
        source=data,
    )
    p1.line(
        "time",
        "local",
        color="#d95f02",
        line_width=1,
        legend_label="Local",
        source=data,
    )
    p1.add_tools(
        HoverTool(
            tooltips=[
                ("Time", "$x"),
                ("Baseline", "@baseline"),
                ("Local", "@local"),
                ("Normalized Error", "@error"),
            ],
            mode="vline",
        )
    )

    # Plot the normalized difference between local and baseline
    p2 = figure(title=title2, x_range=p1.x_range)
    p2.line(
        "time",
        "error",
        color="#7570b3",
        line_width=1,
        legend_label="Normalized Error",
        source=data,
    )
    p2.add_tools(
        HoverTool(tooltips=[("Time", "$x"), ("Normalized Error", "$y")], mode="vline")
    )

    # Shared formatting
    for _plot in (p1, p2):
        _plot.title.align = "center"
        _plot.xaxis.axis_label = xlabel
        _plot.grid.grid_line_alpha = 0.3

    grid = gridplot(
        [[p1, p2]], plot_width=650, plot_height=375, sizing_mode="scale_both"
    )
    script, div = components(grid)
    return script, div


def plot_error(
    baseline_data: list, test_data: list, attributes: List[Tuple[str, str]],
) -> List[Tuple[str, str]]:
    """
    Plots the raw baseline vs test results for each attribute in one column and
    the normalized differences between the two in a second column.

    Parameters
    ----------
    baseline_data : list
        Baseline output included in distribution.
    test_data : list
        User-produced output from OpenFAST.
    attributes : List[Tuple[str, str]]
        List of tuples of attribute names and units.

    Returns
    -------
    List[Tuple[str, str, str]]
        List of tuples of script, div, and attrbiute name for each attribute in
        `attributes`.
    """

    x = test_data[:, 0]
    plots = []
    for i, (name, units) in enumerate(attributes):

        title1 = f"{name} ({units})"
        title2 = "Normalized Difference"
        xlabel = "Time (s)"

        y1 = np.array(baseline_data[:, i], dtype=np.float)
        y2 = np.array(test_data[:, i], dtype=np.float)

        script, div = plot_single_attribute_error(x, y1, y2, xlabel, title1, title2)
        plots.append((script, div, name))

    return plots


def create_plot_body(html_head: str, plots: List[tuple]):
    """
    Updates the HTML header with the required Bokeh scripts and creates a div
    body to be inserted at the end of the HTML file.

    Parameters
    ----------
    html_head : str
        Stringified version of the HTML head.
    plots : List[tuple]
        script, div, attribute tuple for each of the attributes to be included.

    Returns
    -------
    html_head : str
        Updated HTML header.
    div_body : str
        Div objects to be included in the HTML body.
    """

    if not plots:
        return html_head, ""

    script_ix = html_head.rfind("</script>\n") + len("</script>\n")

    div_body = ""
    script_body = "\n"
    for script, div, attribute in plots:
        div = _replace_id_div_string(div, attribute)
        div_body = "".join((div_body, div))

        script = _replace_id_script_string(script, attribute)
        script_body = "\n".join((script_body, script))

    html_head = script_body.join((html_head[:script_ix], html_head[script_ix:]))

    return html_head, div_body


def create_case_summary(  ###### NEED TO ACTUALLY CREATE THE PLOTS
    path: str,
    case: str,
    results: np.ndarray,
    results_max: np.ndarray,
    attributes: List[Tuple[str, str]],
    results_columns: List[str],
    plots: List[Tuple[str, str, str]],
    tolerance: float,
):
    """
    Creates the case summary and exports it to `path`/`case`_summary.html.

    Parameters
    ----------
    path : str
        Path for where to save the html file.
    case : str
        Name of the case.
    results : np.ndarray, shape: [n_attributes, n_norms]
        Norm results by attribute.
    results_max : array-like
        Index for the max value of each norm.
    attributes : List[Tuple[str, str]]
        List of tuples of attribute names and units.
    results_columns : List[str], optional
        List of norms that are being provided, by default ["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"]
    plots : List[Tuple[str, str, str]]
        List of tuples of scipt, div, and attribute name.
    tolerance : float
        Value at which all normed values must be below.
    """

    title = " ".join((case, "Summary"))
    html_head = create_header(title)
    columns = ["Channel", *[r.replace("_", " ").title() for r in results_columns]]
    table_head = create_table_head(columns)

    data = [
        (f'<a href="#{attribute}">{attribute}</a>', *norms)
        for (attribute, _), *norms in zip(attributes, results)
    ]
    table_body = ""
    for i, d in enumerate(data):
        table_body = "\n".join(
            (
                table_body,
                f"{INDENT * 3}<tr>",
                # f'{INDENT * 4}<th scope="row">{i + 1}</th>',
                f"{INDENT * 4}<td>{d[0]}</td>",
            )
        )
        for j, val in enumerate(d[1]):
            if i == results_max[j]:
                _class = ' class="cell-highlight"'
            elif val > tolerance:
                _class = ' class="cell-warning"'
            else:
                _class = ""

            cell = f"{INDENT * 4}<td{_class}>{val:0.4e}</td>"

            table_body = "\n".join((table_body, cell))
        table_body = "\n".join((table_body, f"{INDENT * 3}</tr>"))
    table_body = "\n".join((table_body, f"{INDENT * 2}</table>"))

    html_head, plot_body = create_plot_body(html_head, plots)

    html = "\n".join(
        (
            html_head,
            "",
            "<body>",
            f'{INDENT}<h2 class="text-center">{title}</h2>',
            f'{INDENT}<h4 class="text-center">Maximum values for each norm are <span class="cell-highlight">highlighted</span> and failing norms (norm >= {tolerance}) are <span class="cell-warning">highlighted</span></h4>',
            f'{INDENT}<div class="container">',
            table_head,
            table_body,
            f"{INDENT * 2}<br>",
            f"{INDENT}</div>",
            plot_body,
            # f"{INDENT * 2}</div>",
            # f"{INDENT}</div>",
            "</body>",
            create_tail(),
        )
    )
    with open(os.path.join(path, ".".join((case, "html"))), "w") as f:
        f.write(html)


# Left off here


# def exportResultsSummary(path, results):
#     with open(os.path.join(path, "regression_test_summary.html"), "w") as html:

#         html.write(_htmlHead("Regression Test Summary"))

#         html.write("<body>" + "\n")
#         html.write(
#             '  <h2 class="text-center">{}</h2>'.format("Regression Test Summary") + "\n"
#         )
#         html.write('  <div class="container">' + "\n")

#         # Test Case - Pass/Fail - Max Relative Norm
#         data = [
#             ('<a href="{0}/{0}.html">{0}</a>'.format(r[0]), r[1])
#             for i, r in enumerate(results)
#         ]
#         table = _tableHead(["Test Case", "Pass/Fail"])
#         body = "      <tbody>" + "\n"
#         for i, d in enumerate(data):
#             body += "        <tr>" + "\n"
#             body += '          <th scope="row">{}</th>'.format(i + 1) + "\n"
#             body += "          <td>{0:s}</td>".format(d[0]) + "\n"

#             fmt = "{0:s}"
#             if d[1] == "FAIL":
#                 body += ('          <td class="cell-warning">' + fmt + "</td>").format(
#                     d[1]
#                 ) + "\n"
#             else:
#                 body += ("          <td>" + fmt + "</td>").format(d[1]) + "\n"

#             body += "        </tr>" + "\n"
#         body += "      </tbody>" + "\n"
#         table += body
#         table += "    </table>" + "\n"
#         html.write(table)

#         html.write("    <br>" + "\n")
#         html.write("  </div>" + "\n")
#         html.write("</body>" + "\n")
#         html.write(_htmlTail())
#     html.close()
