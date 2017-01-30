How to use notebooks in IC
==========================

The problem with notebooks
--------------------------

Notebooks are as very useful tool for prototyping and preparing analysis, but
they can create conflicts under a version control system (like
git). Usually this is due to the outputs, which can change from one execution to
another (simply the execution counter of each notebook cell can give
conflicts).

To solve those issues in ICARO we will only put into git notebooks without output. We use
`git filters
<http://pascalbugnion.net/blog/ipython-notebooks-and-git.html>`_ to automatize
this process. This solution requires a filter script (``ipynb_drop_output``) and
some git configuration: the files ``.gitattributes`` and ``.gitconfig`` already
included in the repository and a small change in ``.git/config`` which can be
done with this command:

``git config --add include.path $ICTDIR/.gitconfig``


Another source of trouble when working with notebooks are filenames, since each developer can have different
paths and this will require changes in the code to be able to read them. The
easiest way to avoid this conflicts is by following a naming convention, in our
case **all files must be under $IC_DATA**. Each user have to define ``$IC_DATA``
to some folder and put the files needed in the notebook there.

Workflow
--------

1. First thing you will need to do is to configure your environment,  this can
   be done running the script manage.sh. For example (there are other options) ``source manage.sh  work_in_python_version 3.5``.

2. Define your ``$IC_DATA`` if you don't have it yet.

3. Create/checkout your development branch for the notebook (see ICARO)

4. Once you have it, you can simply run ``jupyter notebook`` and get started
   working on your new notebook. Remember to use ``$IC_DATA`` for your input
   files in case you need them.

5. You can commit your work and push it to you fork as many times you want, but
   keep in mind that the committed version will not have the output as they will
   pass through the stripping filter.

6. In principle that should be enough, if you simply run a notebook without
   changing the code, no changes should be detected by git. But they will
   usually appear if you run ``git status``, for example:

::

 $ git status
 On branch notebooks
 Changes not staged for commit:
   (use "git add <file>..." to update what will be committed)
   (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   ../invisible_cities/cities/diomira.ipynb


If you try to see the changes, you won't be able to see them, ``git diff`` will
not show any output. `This has been reported by other git users
<http://stackoverflow.com/questions/19807979/why-does-git-status-ignore-the-gitattributes-clean-filter>`_. To
solve it, git index has to be updated, the easiest way to do this is running
``git add notebook.ipynb``.


In magit, the notebook won't appear in ``magit-status`` but if you try to change
to another branch, rebase, etc. git will complain:


``GitError! The following untracked working tree files would be overwritten by
checkout:  [Type `$' for details]``

By pressing ``$`` you can see the error:

::

   1 git … checkout master
    error: Your local changes to the following files would be overwritten by checkout:
	invisible_cities/cities/diomira.ipynb
   Please, commit your changes or stash them before you can switch branches.

You can add the file typing ``s`` for the file stage menu and then putting the
notebook giving the conflict. After that git won't detect the changes and you
will be able to change branch or whatever you need.

The Gallery
------------

Since the NBs you are committing to git do not include the output, it is not possible to see plots, fits or results of calculations by simply inspecting them.

On the other hand, since your NB may use input files, any one trying to execute it has to have access to those files.

The *Gallery* is an file and html server, whose roles is:

  1) To allow you to push there the .html versions of your NB (thus full output).
  2) To allow you to push there the files that you are using in your NB so that other users can download them in order to run your script.

ICARO includes three python scripts to facilitate communication with the Gallery. **upload_nb.py** will upload an html copy of your NB to the gallery, **upload_files.py** will upload a list of files from IC_DATA to the gallery and **download_files.py** will download a list of files from the gallery to IC_DATA. 
