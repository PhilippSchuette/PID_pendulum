/*
 * main.js is part of PID_pendulum python module.
 *
 * Author: Daniel Schuette
 * License: GPL-3.0
 * Date: 13/05/2019
 */
/*
 * TODO:
 * - add input fields for pendulum parameters
 * - bootstrap the website layout
 * - add loading animation while data is not fetched from remote origin
 * - save data to buffer and only fetch if absolutely necessary to reduce lag
 * - make pendulum animation
 */
const baseURL = "http://localhost:5000";
const APIRoute = "/api/v1/";
const request = "?alpha=4.0&beta=1.5&mu=0.8&phi0=0.25&phi0_dot=0.25&max_control=3.0&frequency=10&deadband=0.01&set_point=0&precision=5&key=0&t_end=";
const address = baseURL + APIRoute + request;

function convert(data) {
    let res = [];
    for (let i=0; i<data.angles.length; i++) {
        res.push({
            angle: data.angles[i],
            time: data.support_values[i],
        });
    }
    console.log(data);
    console.log(res);
    return res;
}

/* Requires an SVG element to exist. */
function plot(data) {
    const svgWidth = 900;
    const svgHeight= 600;
    const margins = { top: 20, right: 20, bottom: 30, left: 50 };
    const width = svgWidth - margins.left - margins.right;
    const height = svgHeight - margins.top - margins.bottom;

    // clear svg and all child nodes, then recreate svg
    d3.selectAll("svg > *").remove();
    let svg = d3.select("svg")
        .attr("width", svgWidth)
        .attr("height", svgHeight);
    let graph = svg.append("g")
        .attr("transform", "translate("+margins.left+","+margins.top+")");

    let x = d3.scaleLinear().rangeRound([0, width]);
    let y = d3.scaleLinear().rangeRound([height, 0]);
    let line = d3.line()
        .x((d) => x(d.time))
        .y((d) => y(d.angle));
    x.domain(d3.extent(data, (d) => d.time));
    y.domain(d3.extent(data, (d) => d.angle));

    graph.append("g")
        .attr("transform", "translate(0,"+height+")")
        .call(d3.axisBottom(x))
        .select(".domain");

    graph.append("g")
        .call(d3.axisLeft(y))
        .append("text")
        .attr("fill", "#000")
        .attr("font-size", "14px")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", "-3.1em")
        .attr("dx", -height/2+24)
        .attr("text-anchor", "end")
        .text("Angle");

    graph.append("path")
        .datum(data)
        .attr("fill", "none")
        .attr("stroke", "steelblue")
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke-width", 1.5)
        .attr("d", line);
}

window.onload = () => {
    /*
     * TODO: collect data into a buffer and only re-fetch if
     *       the buffer is empty
     * TODO: get all magic numbers from UI
     */
    const animationSpeed = 200; /* one animation per `x' millisecs */
    const numFrames = 100;      /* `x' animations total */

    for (let i=1; i<numFrames; i++) {
        setTimeout(() => {
        fetch(address+i, { method: "get" })
            .then(resp => resp.json())
            .then(json => { plot(convert(json)); })
            .catch(err => {
                console.log(`Error: ${err}`);
                return;
            });
        }, animationSpeed*i);
    }
}
