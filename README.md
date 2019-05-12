# PID pendulum controller
## Overview
This is a short description of the PID controlled pendulum project. The project illustrates several concepts of real-world PID controllers, that theoretical descriptions often disregard. At the moment this includes

- bounded controller output,

- limited controller speed,

- controller measurement imprecision,

- incomplete system information.

Additional features will be added in the future.

## Web API and Demonstration
A demonstration of the PID pendulum controller implementation will be available soon. To run the web app locally, `pip install flask` and run:

```bash
cd api
flask run # now, visit localhost:5000 in a web browser
```

If the webserver is running, API requests can be send to the route `/api/v1`. The following script demonstrates that:

```bash
./req_api.sh # performs a valid API request
./req_api.sh err # performs an invalid API request
```

Currently, the API simply echos the given parameters or returns an error if a required parameter was not provided. Soon, the API will return `json` containing the function values that correspond to the user-provided arguments.

## License
The code in this repository is [GPL-3.0 licensed](./LICENSE.md). `html` and `css` files are not distributed under any license if not specified otherwise.
