"""Provides the error plotting and html generation functionality."""


import os
from typing import List, Tuple

import numpy as np
from bokeh.embed import components
from bokeh.layouts import gridplot
from bokeh.plotting import figure
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
    return div


def _tidy_div_string(div: str) -> str:
    """
    Tidies the stringified div object to have the same class and style.

    Parameters
    ----------
    div : str
        Div object as a string.

    Returns
    -------
    str
        Div object as a string that has been tidied up.
    """

    div_class = ' class="col-sm-12 col-md-6 col-lg-6"'
    style = 'style="margin:10 auto"'

    ix_insert = div.find("></div>")
    div = div_class.join((div[:ix_insert], div[ix_insert:]))
    div = div.replace("<div", " ".join(("<div", style)))
    div = div.join((INDENT * 3, "\n"))
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
    return script


def _tidy_script_string(script: str) -> str:
    """
    Replaces "\n" with a "\n  ".

    Parameters
    ----------
    script : str
        Text of html script.

    Returns
    -------
    str
        Updated text of html script.
    """
    return script.replace("\n", f"\n{INDENT}")


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
            f"{INDENT}<title>{title}</title>'.format()",
            f'{INDENT}<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">',
            f'{INDENT}<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>',
            f'{INDENT}<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>',
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

    header = "\n".join((f"{INDENT * 5}<th>{col}</th>" for col in columns))

    head = "\n".join(
        (
            f'{INDENT * 2}<table class="table table-bordered table-hover table-sm" style="margin: auto; width: 100%; font-size:80%">',
            f"{INDENT * 3}<thead>",
            f"{INDENT * 4}<tr>",
            f"{INDENT * 5}<th>#</th>",
            header,
            f"{INDENT * 4}</tr>",
            f"{INDENT * 3}</thead>",
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

    # Plot the unadjusted comparison of values
    p1 = figure(title=title1)
    p1.line(x, y1, color="#1b9e77", linewidth=3, legend="Baseline")
    p1.line(x, y2, color="#d95f02", linewidth=1, lengend="Local")
    p1.add_tools(HoverTool(tooltips=[("Time", "$x"), ("Value", "$y")]), mode="vline")

    # Plot the normalized difference between local and baseline
    norm = np.absolute(y1 - y2) / y1
    p2 = figure(title=title2)
    p2.line(x, norm, color="#7570b3", linewidth=1, legend="Normalized Error")
    p2.add_tools(
        HoverTool(tooltips=[("Time", "$x"), ("Normalized Error", "$y")], mode="vline")
    )

    # Shared formatting
    for _plot, _title in (p1, p2):
        _plot.title.align = "center"
        _plot.xaxis.axis_label = xlabel
        _plot.grid.grid_line_alpha = 0.3

    grid = gridplot(
        [[p1, p2]], plot_width=650, plot_height=375, sizing_mode="scale_both"
    )
    script, div = components(grid)
    return script, div


def plot_error(
    baseline_data: list,
    baseline_info: list,
    test_data: list,
    test_info: list,
    attributes: list,
) -> List[Tuple[str, str, str]]:
    """
    Plots the raw baseline vs test results for each attribute in one column and
    the normalized differences between the two in a second column.

    Parameters
    ----------
    baseline_data : list
        Baseline output included in distribution.
    baseline_info : list
        Attribute information for `baseline_data`.
    test_data : list
        User-produced output from OpenFAST.
    test_info : list
        Attribute information for `test_data`.
    attributes : list
        list of attribute names.

    Returns
    -------
    List[Tuple[str, str, str]]
        List of tuples of script, div, and attrbiute name for each attribute in
        `attributes`.
    """

    x = test_data[:, 0]
    plots = []
    for attribute in attributes:
        channel = test_info["attribute_names"].index(attribute)

        title1 = "".join((attribute, " (", test_info["attribute_units"][channel], ")"))
        title2 = "Normalized Difference"
        xlabel = "Time (s)"

        y1 = np.array(baseline_data[:, channel], dtype=np.float)
        y2 = np.array(test_data[:, channel], dtype=np.float)

        script, div = plot_single_attribute_error(x, y1, y2, xlabel, title1, title2)
        plots.append((script, div, attribute))

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

    script_ix = html_head.rfind("</script>\n") + len("</script>\n")

    div_body = "\n"
    script_body = "\n"
    for script, div, attribute in plots:
        div = _tidy_div_string(_replace_id_div_string(div, attribute))
        div_body = "\n".join((div_body, div))

        script = _tidy_script_string(_replace_id_script_string(script, attribute))
        script_body = "\n".join((script_body, script))

    html_head = script_body.join((html_head[:script_ix], html_head[script_ix:]))

    return html_head, div_body


def create_case_summary(  ###### NEED TO ACTUALLY CREATE THE PLOTS
    path: str,
    case: str,
    results: np.ndarray,
    results_max: np.ndarray,
    tolerance: float,
    plots: List[str],
    results_columns: List[str] = [
        "max_norm",
        "max_norm_over_range",
        "l2_norm",
        "relative_l2_norm",
    ],
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
        Max of each norm.
    tolerance : float
        Value at which all normed values must be below.
    plots : List[str]
        List of attributes to create plots for.
    results_columns : List[str], optional
        List of norms that are being provided, by default ["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"]
    """

    title = " ".join((case, "Summary"))
    html_head = create_header(title)
    columns = ["Channel", *[r.replace("_", " ").title() for r in results_columns]]
    table_head = create_table_head(columns)

    data = [
        (f'<a href="#{attribute}">{attribute}</a>', *norms)
        for attribute, *norms in results
    ]

    table_body = "".join((INDENT * 3, "<tbody>"))
    for i, d in enumerate(data):
        table_body = "\n".join(
            (
                table_body,
                f"{INDENT * 4}<tr>",
                f'{INDENT * 5}<th scope="row">{i + 1}</th>',
                f"{INDENT * 5}<td>{d[0]}</td>",
            )
        )
        for j, val in enumerate(d[1]):
            if val == results_max[j]:
                _class = ' class="cell-warning"'
            elif val > tolerance:
                _class = ' class="cell-highlight"'
            else:
                _class = ""

            cell = f"{INDENT * 5}<td{_class}>{val:0.4e}</td>"

            table_body = "\n".join((table_body, cell))
        table_body = "\n".join((table_body, f"{INDENT * 4}</tr>"))
    table_body = "\n".join((f"{INDENT * 3}</tbody>", f"{INDENT * 2}</table>"))

    html_head, plot_body = create_plot_body(html_head, plots)
    if not plots:
        plot_body = ""

    html = "\n".join(
        (
            html_head,
            "",
            "<body>",
            f'{INDENT}<h2 class="text-center">{title}</h2>',
            f'{INDENT}<h4 class="text-center">Maximum values for each norm are <span class="cell-warning">highlighted</span> and failing norms (norm >= {tolerance}) are <span class="cell-highlight">highlighted</span></h2>',
            f'{INDENT}<div class="container"',
            table_head,
            table_body,
            f"{INDENT * 2}<br>",
            f"{INDENT}</div>",
            plot_body,
            f"{INDENT * 2}</div>",
            f"{INDENT}</div>",
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
