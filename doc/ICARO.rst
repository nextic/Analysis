ICARO
=====

ICARO stands for Invisible Cities Analysis RepO (yes, we strech the acronym a bit, and yes, we like to fly near the sun too). ICARO is a *client* of IC, the *Invisible Cities* python-based software package for NEXT. The *input* for ICARO analysis are files (on hdf5 format), containing high-level objects produced by IC, such as PMAPS, tracks, Paolina-tracks, and so on. The available tools are those provided by the IC ecosystem (Numpy, SciPy, matplotlib, pytables, pandas...) plus all the utilities provided by IC itself (access to data base, core functions, etc).

Like IC, the programming language for ICARO is python (and cython). While not compulsory, the analysis in ICARO are expected to use heavily the newly available technology of Jupyter notebooks to draft, develop and document analysis. However, it is also expected that production-ready analysis are eventually expressed as a configurable python script. Neither IC, nor ICARO provide a rigid "framework" in which the an analyst can sandwich her "user functions" (we believe in freedom), but rather offers tools and guidance for the analysis to write her own code.

A convenient way to express an analysis script is to prepare it as a "state machine". A SM can be viewed as an automaton programmed to run a task by previously setting a set of switches. Before running the task, the SM checks that all switches are "on" and therefore its configuration ready. Typically, a SM is set by defining input and output files, printing an plotting options, and giving values to constants needed in order to perform the automaton tasks. For two examples of SM, see the cities *diomira* and *irene* in IC.

Before writing an analysis script, however, the analyist can draft her work in a Jupyter Notebook (NB). NBs are very dynamical by nature and allow mixing code, rich text and plots, thus facilitating enormously the task of drafting (and later documenting) an analysis. See the section on NBs in this docs.

 The recommended organization for ICARO is one folder per analysis (e.g, NA-22, Co-56, Krypton, bb2nu...). Additional folders can be created inside the specific analysis folders if needed. In addition to analysis folders, several other directories are created. Specifically:

1) Directory **/cities** containing notebooks (NBs) which document specific *cities* belonging to *IC*. The notebooks are included in ICARO as an example of documentation NBs.
2) Directory **/calib** containing NBs and scripts for calibration algorithms.
3) Directory **/utils** containing NBs documenting utilities for analysis.

Notice that functions and classes whose scope is not specific to a given analysis (e.g, a fitting function), will normally be exported to the server repository (IC).

Also, notice that the main (central or upstream) ICARO repository is expected to contain only mature NBs and scripts. Notebooks tend to proliferate in the development phase of an analysis. This is perfectly OK, provided that the analyst keeps the little monsters in her own branch and only commits upstream those mature enough to describe/document her analysis (or analysis step).

ICARO Workflow
===============

To illustrate the work flow in ICARO, consider developer Leia, working on her brand new analysis *dark_vader*

1) Leia forks the main (upstream) ICARO repository, which is expected to contain always a snapshot of a "mature state" of analysis. This does not mean that Leia and the other analysts (Hans, Luke, etc.) wait until their analysis are fully finished to commit a NB or script upstream. Instead, the idea is to commit NBs that describe reasonably well understood steps/stages/methods in the analysis. This, often, requires working for a while with a draft and fast-evolving NB (or collection of NBs). Those are to stay in Leia (Luke, Hans...) forks.

2) After forking ICARO from upstream, Leia creates a branch *devel* sprouting from *master*.

3) She checkouts *devel*.

4) In *devel* Leia creates a directory **/dark_vader** on her own fork and starts to work there. She commits notebooks happily as she progresses through the difficult problem of analyzing the data needed to destroy the Death Star.

5) At a given point, Leia has a first version of her analysis mature enough to be shared with others. Such analysis can be described, say, in the notebook **/dark_vader/dv.ipynb**. In addition, she has a notebook describing utilities
(e.g, phaser_handling) which she can place in **/utils/phaser_handling.ipynb**, since this can be useful for other analysis. Furthermore, she has written a few nice fitting functions which are heavily use by **dv.ipynb**. Finally, she has a large number of drafts (**dv_1.ipynb**, **dv_2.ipynb**...) hanging around on her *devel* branch in her own fork.

6) In order to ship her work to the rest of the team, Leia needs to make a Pull Request (PR) upstream. The steps are:
    a) She makes sure that her *main* is up-to-date with *upstream/main* (resetting, if needed her main to upstream/main).
    b) She creates a branch *release* sprouting from *main* (if she has already a running *release* branch she resets it to *main*).
    c) She checks out into *release* the notebook(s) and/or script(s) she wants to commit upstream. (using magit-checkout-file), e.g, **/dark_vader/dv.ipynb**, **/utils/phaser_handler.ipynb**.
    d) She also needs to commit **fitting_functions.py** to /core in IC (more on that later).
    d) She makes a pull request.

Notice that the workflow in ICARO is different from the workflow in IC, where the developers are normally continuously rebasing their *master* to *upstream/master* and sprouting short-level branches from there. In ICARO, instead, the analysis keeps her development in the long-lived branch *devel* (she can sprout as many branches as she likes from there). From *devel* she selects the (few) notebooks and scripts that need to go upstream, through the branch *release*. Thus, in ICARO, *main* and *release* (which can be a long running branch following *master* or a short lived branch sprouting from *main* created just for the specific purpose of delivering goods upstream) follow the state of the upstream repository, while *devel* evolves with the work of the analyst.

To commit **fitting_functions.py** to IC, Leia can proceed in two different ways.

 a) If she is an IC contributor and has her own fork of IC, she makes an ordinary PR to IC including her **fitting_functions.py**

 b) If she is not an IC contributor, she includes her **fitting_functions.py** in her PR to ICARO, in this case in the directory **/ic_deliver/core/fitting_functions.py**. Notice that one also expects a **fitting_functions_test.py**. Once her PR is approved, the dictator of ICARO will move the code to IC.

ICARO Policy
=============

 ICARO is running by a dictator (Daedalus), who is in charge of coordinating the reviewing and approval of all the code committed to ICARO. When a PR is submitted to ICARO, Daedalus will assign the revision of the PR to one of the active contributors in ICARO (including himself), and/or and external reviewer.

 A PR review on ICARO *always* imply that the reviewer can run the PR code, with a sample data. To facilitate this, the PR should include a NB describing and exercising the analysis. The data used by this PR must be uploaded to the *gallery* (see separate note) before the PR is submitted. The data in the NB is always pointed to by the environment variable $IC_DATA. In order to run the NB, the reviewer will access the source code through git, download the needed files from the gallery to her own $IC_DATA area and execute the NB. She can then comment on the analysis and eventually approve it.

 As an example, assume that Leia submits a PR upstream containing her notebook **dv.ipynb**. The NB reads a file
 **ds_tag.h5** which must be in Leia's $IC_DATA area. Notice the "tag" in the data file which stands for a tag, hash, or label (e.g, a stream describing the time or a randomly generated hash) that makes the name unique. The steps are:

  1) Leia uploads **ds_tag.h5** to the gallery.
  2. Leia submits a PR wich includes **dv.ipynb**.
  3. Daedalus assigns Jyn as a reviewer.
  4. Jyn pulls Leia fork branch and gets **dv.ipynb**. She also download the file **ds_tag.h5** from the gallery to her own $IC_DATA area.
  5. Jyn runs the notebook, and comments on results. She can eventually propose improvements to Leia, and both of them can work on the same NB through mutual PR across their own forks.
  6. Eventually they converge, Jyn write her positive report and Daedalus approves the PR. **dv.ipynb** is now part of ICARO.

Notice that this model favour *pair programming*, through the revision process. The idea is that all the analysis in ICARO are understood *quantitatively* by more than one analyst, in order to minimize errors and maximize exchange of ideas.

Happy flying! May the force be with you! 
