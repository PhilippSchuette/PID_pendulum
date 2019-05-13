# PID pendulum controller
## Overview
This is a short description of the PID controlled pendulum project. The project illustrates several concepts of real-world PID controllers, that theoretical descriptions often disregard. At the moment this includes

- bounded controller output,

- limited controller speed,

- controller measurement imprecision,

- incomplete system information.

Additional features will be added in the future.

## Build Sources and Documentation
Make sure to install the following dependencies with Python 3 (might require `pip3` on your system):

- `pip install flask` for the web app
- `pip install flake8` for linting
- `pip install sphinx` for all documentation

You might also need:
- `make` for build from source
- `LaTeX` for creating PDF documentation

Then, run:

```bash
make all # install the PID_pendulum package
make # show a list of documentation formats, e.g. `make latexpdf` creates PDF docs
```

You can now import the pendulum in your code as follows:

```python
from PID_pendulum.PID_control import Pendulum, PIDControl
# see documentation for pendulum and controller methods
```

## Web API and Demonstration
An online demonstration of the PID pendulum controller implementation will be available soon. To run the web app locally though, run:

```bash
cd api
flask run # now, visit localhost:5000 in a web browser
```

If the webserver is running, API requests can be send to the route `/api/v1`. The following script demonstrates that:

```bash
cd utils
./req_api.sh # performs an API request with valid parameters
./req_api.sh err # performs an API request with invalid parameters
```

Currently, the API simply echos the given parameters or returns an error if a required parameter was not provided. Soon, the API will return `json` containing the function values that correspond to the user-provided arguments.

## License
The code in this repository is [GPL-3.0 licensed](./LICENSE.md). `html` and `css` files are not distributed under any license if not stated otherwise in the files themselves.
