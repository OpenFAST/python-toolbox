"""Provides the error plotting and html generation functionality."""


import os
import sys

import bokeh
import numpy as np
from bokeh.embed import components
from bokeh.layouts import gridplot
from bokeh.models.tools import HoverTool
from bokeh.models.widgets import RangeSlider
from bokeh.plotting import figure

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


def _tidy_script_string(script):
    return script.replace("\n", "".join(("\n", INDENT)))


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
            "".join((INDENT, f"<title>{title}</title>'.format()")),
            "".join(
                (
                    INDENT,
                    '<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">',
                )
            ),
            "".join(
                (
                    INDENT,
                    '<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>',
                )
            ),
            "".join(
                (
                    INDENT,
                    '<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>',
                )
            ),
            "".join((INDENT, '<style media="screen" type="text/css">')),
            "".join((INDENT * 2, ".cell-warning {")),
            "".join((INDENT * 3, "background-color: #FF6666;")),
            "".join((INDENT * 2, "}")),
            "".join((INDENT * 2, ".cell-highlight {")),
            "".join((INDENT * 3, "background-color: #E5E589;")),
            "".join((INDENT * 2, "}")),
            "".join((INDENT, "</style>")),
            "</head>",
        )
    )
    return html


def create_tail():
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

    header = "\n".join(
        ("".join((INDENT * 5, f"<th>{column}</th>")) for column in columns)
    )

    head = "\n".join(
        (
            "".join(
                (
                    INDENT * 2,
                    '<table class="table table-bordered table-hover table-sm" style="margin: auto; width: 100%; font-size:80%">',
                )
            ),
            "".join((INDENT * 3, "<thead>")),
            "".join((INDENT * 4, "<tr>")),
            "".join((INDENT * 5, "<th>#</th>")),
            header,
            "".join((INDENT * 4, "</tr>")),
            "".join((INDENT * 3, "</thead>")),
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
) -> list:

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
    for i, (script, div, attribute) in enumerate(plots):
        div = _tidy_div_string(_replace_id_div(div, attribute))
        div_body = "\n".join((div_body, div))

        script = _tidy_script_string(_replace_id_script_string(script, attribute))
        script_body = "\n".join((script_body, script))

    html_head = script_body.join((html_head[:script_ix], html_head[script_ix:]))

    return html_head, div_body


def create_case_summary(path, case, results, results_max, tolerance, plots):

    title = " ".join((case, "Summary"))
    html_head = create_header(title)
    html_head, div_body = create_plot_body(html_head, plots)
    columns = [
        "Channel",
        "Max Norm",
        "Relative Max Norm",
        "L2 Norm",
        "Relative L2 Norm",
    ]
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
                "".join((INDENT * 4, "<tr>")),
                "".join((INDENT * 5, f'<th scope="row">{i + 1}</th>')),
                "".join((INDENT * 5, f"<td>{d[0]}</td>")),
            )
        )
        for j, val in enumerate(d[1]):
            if val == results_max[j]:
                _class = ' class="cell-warning"'
            elif val > tolerance:
                _class = ' class="cell-highlight"'
            else:
                _class = ""

            cell = "".join((INDENT * 5, f"<td{_class}>{val:0.4e}</td>"))

            table_body = "\n".join((table_body, cell))
        table_body = "\n".join(("".join(INDENT * 4, "</tr>")))
    table_body = "\n".join(
        ("".join((INDENT * 3, "</tbody>")), "".join((INDENT * 2, "</table>")))
    )

    if len(plots) == 0:
        plot_body = ""
    else:
        html_head, plot_body = create_plot_body(html_head, plots)

    html = "\n".join(
        (
            html_head,
            "",
            "<body>",
            "".join((INDENT, f'<h2 class="text-center">{title}</h2>')),
            "".join(
                (
                    INDENT,
                    f'<h4 class="text-center">Maximum values for each norm are <span class="cell-warning">highlighted</span> and failing norms (norm >= {tolerance}) are <span class="cell-highlight">highlighted</span></h2>',
                )
            ),
            "".join((INDENT, '<div class="container"')),
            table_body,
            "".join((INDENT * 2, "<br>")),
            "".join((INDENT, "</div>")),
            plot_body,
            "".join((INDENT * 2, "</div>")),
            "".join((INDENT, "</div>")),
            "</body>",
            create_tail(),
        )
    )
    with open(os.path.join(path, ".".join((case, "html"))), "w") as f:
        f.write(html)


# Left off here


def exportResultsSummary(path, results):
    with open(os.path.join(path, "regression_test_summary.html"), "w") as html:

        html.write(_htmlHead("Regression Test Summary"))

        html.write("<body>" + "\n")
        html.write(
            '  <h2 class="text-center">{}</h2>'.format("Regression Test Summary") + "\n"
        )
        html.write('  <div class="container">' + "\n")

        # Test Case - Pass/Fail - Max Relative Norm
        data = [
            ('<a href="{0}/{0}.html">{0}</a>'.format(r[0]), r[1])
            for i, r in enumerate(results)
        ]
        table = _tableHead(["Test Case", "Pass/Fail"])
        body = "      <tbody>" + "\n"
        for i, d in enumerate(data):
            body += "        <tr>" + "\n"
            body += '          <th scope="row">{}</th>'.format(i + 1) + "\n"
            body += "          <td>{0:s}</td>".format(d[0]) + "\n"

            fmt = "{0:s}"
            if d[1] == "FAIL":
                body += ('          <td class="cell-warning">' + fmt + "</td>").format(
                    d[1]
                ) + "\n"
            else:
                body += ("          <td>" + fmt + "</td>").format(d[1]) + "\n"

            body += "        </tr>" + "\n"
        body += "      </tbody>" + "\n"
        table += body
        table += "    </table>" + "\n"
        html.write(table)

        html.write("    <br>" + "\n")
        html.write("  </div>" + "\n")
        html.write("</body>" + "\n")
        html.write(_htmlTail())
    html.close()
