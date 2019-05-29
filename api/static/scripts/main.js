/*
 * main.js is part of PID_pendulum python module.
 *
 * Author: Daniel Schuette
 * License: GPL-3.0
 * Date: 13/05/2019
 */
/*
 * TODO:
 * - move from d3 to chart.js
 * - add input fields for pendulum parameters
 * - bootstrap the website layout
 * - add loading animation while data is not fetched from remote origin
 * - save data to buffer and only fetch if absolutely necessary to reduce lag
 * - make pendulum animation
 */
const baseUrl = "http://localhost:5000";
const APIRoute = "/api/v1/";
const request = "?alpha=4.0&beta=1.5&mu=0.8&phi0=0.25&phi0_dot=0.25&max_control=3.0&frequency=10&deadband=0.01&set_point=0&precision=5&key=0&t_end=";
const requestUrl = baseUrl + APIRoute + request;

window.onload = () => {
    /*
     * TODO: collect data into a buffer and only re-fetch if
     *       the buffer is empty
     * TODO: get all magic numbers from UI
     */
    const animationSpeed = 200; /* one animation per `x' millisecs */
    const numFrames = 30;       /* `x' animations total */

    for (let i=1; i<numFrames; i++)
        setTimeout(getDataAndApply(requestUrl+i, lineGraph), animationSpeed*i);
};

function lineGraph(data) {
    console.log(data);
}

function getDataAndApply(url, callback) {
    fetch(url, { method: "get" })
        .then(resp => resp.json())
        .then(json => { callback(json); })
        .catch(err => {
            console.log(`Error: ${err}`);
            return;
        });
}
