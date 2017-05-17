
# The "library"

* **models/** — abstract optimizers
* **methods.py** — models, based on optimization methods implementations
* **datetime_functions.py** — auxiliary procedures to work with time format
    given in the dataset

## AbstractBaroyanRvachevOptimizer
To instantiate an object, initializer requires the following arguments:
* *data* = incidence collection (such as list)
* *population_quantity* = number of inhabitants in considered period of time
* *params* = instance of `baroyan_rvachev.BaroyanParams`'s subclass

After that `optimizer.fit_one_outbreak()` void method should be fired to start
computations.

Finally, `optimizer.get_results()` method returns dict(str: object),
which contains all the parameters which provide local optimum *R^2* value.
Set of fields is often changing, but cares about backward-compatibility.

Full example is provided in
[fit_baroyan_rvachev.py](../experiments/fit_baroyan_rvachev.py)
experiment listing.

## AbstractSEIROptimizer
To instantiate an object, initializer requires the following arguments:
* *data* = incidence collection (such as list)
* *initial_data* = initial incidence data before cutting up to `sample_size`.
    Currently disables (equal to `data`), but might be useful in the future.
* *dates* dates (str YYYYMMDD) collection, corresponding to `initial_data`
* *population_quantity* = number of inhabitants in considered period of time
* *params* = instance of `baroyan_rvachev.BaroyanParams`'s subclass

After that `optimizer.fit_one_outbreak()` void method should be fired to start
computations.

Finally, `optimizer.get_results()` method returns dict(str: object),
which contains all the parameters which provide local optimum *R^2* value.
Set of fields is often changing, but cares about backward-compatibility.

Full example is provided in
[fit_seir.py](../experiments/fit_seir.py)
experiment listing.

## Implementations
Next approach is applicable to both optimizers.
Abstract optimizers don't have to be very flexible (at the current stage),
so implementation is obliged to override only one method `optimize`
like in the following example:

```python
class SEIRMyMethodOptimizer(AbstractSEIROptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = ... # minimize the function somehow,
                     # using initial minimize_params
                     # in minimize_params_range bounds
        # resulting fit value, final bunch of optimal values
        return result.fun, tuple(result.x)
```

Few implementations are provided in [methods.py](methods.py).
