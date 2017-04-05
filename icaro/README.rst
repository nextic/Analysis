ICARO
=====

This is the top level of the ICARO repository. Different analysis belong here. The recommended organization is one folder per analysis (e.g, NA-22, Co-56, Krypton, bb2nu...). Additional folders can be created inside the specific analysis folders if needed. In addition to analysis folders, several other directories are created. Specifically:

1) Directory **/cities** containing notebooks (NBs) which document specific **cities** belonging to the reconstruction repository **IC**.
2) Directory **/calib** containing NBs and scripts for calibration algorihtms.
3) Directory **utils** containing NBs documenting utilities for analysis.

Notice that functions and classes whose scope is not specific to a given analysis (e.g, a fitting function), will normally be placed in modules (eg., **fitting_functions.py**) which will be located in the core repository in IC.

Also, notice that the main (central) ICARO repository is expected to contain ony mature NBs and scripts. By its own nature, NBs tend to proliferate in the development phase of an analysis. This is perfectly OK, provided that the developer keeps the NBs in her own branch and only commits to the main repository those mature enough to describe/document her analysis (or analysis step).

To illustrate the work flow in ICARO, consider developer Leia, working on her brand new analysis dark_vader.

1) Leia forks the main ICARO repository, which is expected to contain always a snapshot of a "mature state" of analysis. This does not mean that Leia and the other developers (Hans, Luke, etc.) wait until their analysis are fully finished to commit a NB or script to central repository. Instead, the idea is to commit NBs that describe reasonably well understood steps/stages/methods in the analysis. This, often, requires working for a while with a draft and fast-evolving NB (or collection of NBs). Those are to stay in Leia (Luke, Hans...) branches.

2) She created immediately a branch *devel* and a branch *release*, both sprouting from *master*.

3) She checkout *devel*.

4) In *devel* Leia creates a directory **/dark_vader** on her own fork and starts to work there. She commits NBs happily as she progresses through the difficult problem of analysis the data needed to destroy the Death Star.

5) At a given point, she has a first version of her analysis, mature enough to be shared with others. Such analysis can be described, say, in the notebook **/dark_vader/dv.ipynb**. In addition, she has a notebook describing utilities
(e.g, phaser_handling) which she can place in **/utils/phaser_handling.ipynb**, since this can be useful for other analysis. Furthermore, she has written a few nice fitting functions which are heavily use by **dv.ipynb**. Finally, she has a large number of drafts (dv_1.ipynb, dv_2.ipynb...) hanging around on her devel branch.

6) In order to ship her work to the rest of the team, Leia needs to make a Pull Request (PR). She ten checks out branch *release* and then checks out into that branch the files that she wants to commit, e.g, **/dark_vader/dv.ipynb**,
**/utils/phaser_handler.ipynb**. She also needs to commit fitting_functions.py to /core in IC. She can do this in two different ways:

 a) If she is an IC developer and has her own fork of IC, she makes an ordinary PR to IC including her fitting_functions.py

 b) If she is not an IC developer, she includes her fitting_functions.py in her PR to ICARO, in this case in the directory *ic_deliver/core/fitting_functions.py*. Once her PR is approved, the dictator of ICARO will move her fitting_functions.py to IC.

Notice that this two-branch scheme, allows to keep a very tidy *deliver* branch and a busy-and-messy *devel* branch. Of course, nothing prevents Leia from branching as much as she wants from *devel* itself.

As the central repository grows from commits of Hand, Luke and the others, Leia resets her own *master* (and *deliver*) branches to the central (upstream) master, and rebases her *devel* branch to upstream/master.
