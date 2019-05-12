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

## License
The code in this repository is [GPL-3.0 licensed](./LICENSE.md). `html` and `css` files are not distributed under any license if not specified otherwise.
