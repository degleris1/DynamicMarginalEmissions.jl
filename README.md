# DynamicMarginalEmissions.jl :zap:

The repo implements the dynamic marginal emissions rates calculations derived in:

> **Lucas Fuentes Valenzuela, [Anthony Degleris](https://degleris1.github.io/), [Abbas El Gamal](https://isl.stanford.edu/~abbas/), [Marco Pavone](https://web.stanford.edu/~pavone/), [Ram Rajagopal](https://ramr.su.domains/)**. 
> <br> [Dynamic locational marginal emissions via implicit differentiation](
https://arxiv.org/abs/2302.14282). 
> <br> *IEEE Transactions on Power Systems*. 2023 Feb 22.

The package solves standard dynamic economic dispatch problems used to dispatch electricity systems.
It supports a full electricity network model with linearized (DC) power flow constraints, as well as batteries and ramping constraints.

After computing the optimal dispatch, the package can compute *the derivative of total emissions with respect to changes electricity demand at a given node and time*, known as the **(dynamic) locational marginal emissions rate (LME)**.
The LMEs are **dynamic** :zap: because they not only calculate how changes in demand will affect emissions at the current moment, but also how they will affect emissions later in the day.



## Installation

First install [Julia](https://julialang.org/downloads/). 
Clone this repo by running the following commands:

```
git clone https://github.com/degleris1/DynamicMarginalEmissions.jl.git
```

Finally, navigate to `DynamicMarginalEmissions.jl`, launch `julia`, and run

```
] add .
```




## Usage 

See `notebooks/demo.jl` for the basic functionality of the repo.




## Reproducing paper results

All results generated for our article are in `experiments/rodm` 
and `experiments/rodm`.




## Citation

If you use our work in you research, please cite the following reference.

```
@article{valenzuela2023dynamic,
  title={Dynamic locational marginal emissions via implicit differentiation},
  author={Valenzuela, Lucas Fuentes and Degleris, Anthony and El Gamal, Abbas and Pavone, Marco and Rajagopal, Ram},
  journal={IEEE Transactions on Power Systems},
  year={2023},
  publisher={IEEE}
}
```




## TODO

- [ ] Remove unused exports
- [ ] Add advanced functionality