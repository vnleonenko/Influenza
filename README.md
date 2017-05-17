## Synopsis

TODO at the top of the file there should be a short introduction and/ or overview that explains **what** the project is. this description should match descriptions added for package managers (gemspec, package.json, etc.)

## Code example

TODO show what the library does as concisely as possible, developers should be able to figure out **how** your project solves their problem by looking at the code example. make sure the api you are showing off is obvious, and that your code is short and concise.

## Motivation

TODO a short description of the motivation behind the creation and maintenance of the project. this should explain **why** the project exists.

## Installation

We recommend you to use a virtual environment for this project.
Read more [here][venv-python]

```sh
$ git clone https://github.com/vnleonenko/Influenza
$ cd Influenza
$ python -m venv .venv
$ source .venv/bin/activate  # for Windows ".venv/scripts/activate"
(.venv) $ pip install -r requirements.txt
```

If you have troubles with requirements installation, ensure that you don't
have any spaces in path to your python interpreter. If so, create it in
another place.
If you still have problems, follow installation process for
[linux][venv-linux], [windows][venv-windows] or [os x][venv-osx].

## Project overview

* [**BaroyanRvachev/**](BaroyanRvachev/) — legacy code (deprecated in Nov 2016),
    Baroyan-Rvachev and SEIR models implementations.
    `BaroyanRvachev_v6.py` script still working fine
* [**FluGraphs_v2/**](FluGraphs_v2/) — mostly old code,
    useful visualisation tools for RJNAMM and incidence interpolation
    error computation
* [**core/**](core/) — library, containing abstract Baroyan-Rvachev and
    SEIR models (in **models/** module), and their implementations in
    **methods.py**
* [**data/**](data/) — *incidence* and *population* datasets for
    considered cities - currently {msk, spb, nsk}, provided by
    [Research Institute of Influenza][InfluenzaInstitute]
* [**experiments/**](experiments/) — directory, containing all the experiments and results
    for some of them (included because of huge computational time)

## Tests

TODO describe and show how to run the tests with code examples.
fixme add usage to size_comparison.py how to generate *.png

## Contributors

All the implemented experiments are located in **experiments** directory
and obey the following conventions:

* Paths to results are built using `common.RESULTS_PATH`
* To get population dict for specified city
    `common.get_population(city_mark)` procedure is used
* To obtain a list of filenames with incidence data
    `common.get_incidence_filenames(city_mark)` procedure is used


## License

This project is licensed under the terms of the [GNU GPLv3](LICENSE) license.



[//]: # (these are reference links used in the body of this note and get stripped out when the markdown processor does its job. there is no need to format nicely because it shouldn't be seen. thanks so - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [venv-python]: <https://docs.python.org/3/library/venv.html>
   [venv-linux]: <http://docs.python-guide.org/en/latest/dev/virtualenvs/>
   [venv-windows]: <https://zignar.net/2012/06/17/install-python-on-windows/>
   [venv-osx]: <http://www.marinamele.com/2014/07/install-python3-on-mac-os-x-and-use-virtualenv-and-virtualenvwrapper.html>
   [InfluenzaInstitute]: <http://www.influenza.spb.ru/en/>
   [john gruber]: <http://daringfireball.net>
   [@thomasfuchs]: <http://twitter.com/thomasfuchs>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [ace editor]: <http://ace.ajax.org>