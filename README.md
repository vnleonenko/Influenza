## synopsis

todo at the top of the file there should be a short introduction and/ or overview that explains **what** the project is. this description should match descriptions added for package managers (gemspec, package.json, etc.)

## code example

todo show what the library does as concisely as possible, developers should be able to figure out **how** your project solves their problem by looking at the code example. make sure the api you are showing off is obvious, and that your code is short and concise.

## motivation

todo a short description of the motivation behind the creation and maintenance of the project. this should explain **why** the project exists.

## installation

we recommend you to use virtual environment for this project. read more [here][venv-python]

```sh
$ git clone https://github.com/vnleonenko/influenza
$ cd influenza
$ python -m venv .venv
$ source .venv/bin/activate  # for windows ".venv/scripts/activate"
(.venv) $ pip install -r requirements.txt
```

if you have troubles with requirements installation, ensure that you don't have any spaces in path to your python interpreter. if so, create it in another place.
if you still have problems, follow installation process for [linux][venv-linux], [windows][venv-windows], [os x][venv-osx]

todo provide code examples and explanations of how to get the project.

## project overview

* [**BaroyanRvachev/**](BaroyanRvachev/) — legacy code (deprecated nov 2016),
    baroyan-rvachev and seir models implementations.
    `BaroyanRvachev_v6.py` script still working fine
* [**FluGraphs_v2/**](FluGraphs_v2/) — legacy code (deprecated dec 2016),
    useful visualisation tools for rjnamm and incidence interpolation
    error computation
* **core/** — library, containing abstract baroyan-rvachev and seir models
    (under **models/**), and their implementations in **methods.py**
* **data/** — *incidence* and *population* datasets for
    considered cities - currently msk, spb, nsk, provided by
    [Research Institute of Influenza][InfluenzaInstitute]
* [**experiments/**](experiments/README.md) — directory, containing all the experiments and results
    for some of them (included because of huge computational time)


seirbaroyan/core/ -- корень "библиотеки"
seirbaroyan/core/models/ -- абстрактные оптимизаторы
seirbaroyan/core//methods.py -- реализации методов

seirbaroyan/core/old_code -- старый код, который вы мне отправляли, практически не тронутый. baroyanrvachev_v6.py работает

flugraphs_v2/ -- тоже старый код, незначительные изменения
flugraphs_v2/compile_data_file_v2.py -- вычисление и рисования отклонения интерполированной недельной (суммы дневных( заболеваемости и фактической недельной заболеваемости

## tests

todo describe and show how to run the tests with code examples.
fixme add usage to size_comparison.py how to generate *.png

## contributors

all the implemented experiments are located in **experiments** directory
and obey following conventions for 'low coupling' reasons:

* path to results built using `common.results_path`
* to get population dict for specified city
    `common.get_population(city_mark)` procedure is used
* to obtain a list of filenames with incidence data
    `common.get_incidence_filenames(city_mark)` procedure is used


## license

todo a short snippet describing the license (mit, apache, etc.)



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