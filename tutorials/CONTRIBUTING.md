# Contributing to `dingo`

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to dingo, 
which are hosted in the [GeomScale Organization](https://github.com/GeomScale) on GitHub. 
These are mostly guidelines, not rules. 
Use your best judgment, and feel free to propose changes to this document in a pull request.

## Table of Contents

  * [Prerequisites (how to start)](#prerequisites-how-to-start)
  * [Testing the development branch of `dingo` (get the tools ready)](#testing-the-development-branch-of-dingo-get-the-tools-ready)
  * [Fork  `dingo` repository (this is your repo now!)](#fork-dingo-repository-this-is-your-repo-now)
    + [Verify if your fork works (optional)](#verify-if-your-fork-works-optional)
  * [Working with `dingo` (get ready to contribute)](#working-with-dingo-get-ready-to-contribute)
    + [GitFlow workflow](#gitflow-workflow)
    + [Create new branch for your work](#create-new-branch-for-your-work)
    + [Verify your new branch (optional)](#verify-your-new-branch-optional)
  * [Modify the branch (implement, implement, implement)](#modify-the-branch-implement-implement-implement)
    + [Tests](#tests)
    + [Push](#push)
  * [Pull request (the joy of sharing)](#pull-request-the-joy-of-sharing)
  * [Review (ok this is not an exam)](#review-ok-this-is-not-an-exam)
  
## Prerequisites (how to start)

* git (see [Getting Started with Git](https://help.github.com/en/github/using-git/getting-started-with-git-and-github))
* a compiler to run tests - gcc, clang, etc.
* configured GitHub account
 
Other helpful links:

* http://git-scm.com/documentation
* https://help.github.com/articles/set-up-git
* https://opensource.com/article/18/1/step-step-guide-git

## Testing the development branch of dingo (get the tools ready)

Clone the repository, 

    git clone git@github.com:geomscale/dingo.git dingo
    cd dingo
    git branch -vv

the last command should tell you that you are in `develop` branch.

Now you need to get the `volesti` sumbodule that `dingo` makes use of. 

To do so, you need to run 

    git submodule update --init

Now get the `boost` library:

    wget -O boost_1_76_0.tar.bz2 https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2 
    tar xjf boost_1_76_0.tar.bz2
    rm boost_1_76_0.tar.bz2

And now you are ready to compile `dingo`

    python setup.py install --user 

Once the last command is completed, you may check if everythings is fine by running some `dingo` tests

    python tests/unit_tests.py


If everything is ok, you will see something like this:

[![asciicast](https://asciinema.org/a/3IwNykajlDGEndX2rUtc0D2Ag.svg)](https://asciinema.org/a/3IwNykajlDGEndX2rUtc0D2Ag)

## Fork `dingo` repository (this is your repo now!)

You can't work directly in the original `dingo` repository, therefore you should create your fork of this library. 
This way you can modify the code and when the job is done send a pull request to merge your changes with the original 
repository.

![fork](https://raw.githubusercontent.com/hariszaf/dingo/dingo_tutorial/tutorials/figs/fork.png)

1. login on `GitHub`
2. go to [dingo repository](https://github.com/GeomScale/dingo)
3. click the 'Fork' button
4. choose your profile
5. wait
6. ready to contribute!

More info: [Forking Projects](https://guides.github.com/activities/forking/)

### Verify if your fork works (optional)

Go out of `dingo` directory

    cd ..

clone your repository and checkout develop branch

    git clone git@github.com:hariszaf/dingo.git dingo_fork
    cd dingo_fork
    git checkout develop
    git branch -vv
    git pull

In this case `hariszaf` was the user's name who had forked `dingo`. Make sure you replace this with your own. 

To see the so far commits, simply run:

    git log
    gitk

For now you should see exactly the same commits as in `dingo` repository.

## Working with `dingo` (get ready to contribute)

### GitFlow workflow

`dingo` is using the [GitFlow](http://nvie.com/posts/a-successful-git-branching-model/) workflow. 
It's because it is very well suited to collaboration and scaling the development team. 
Each repository using this model should contain two main branches:

* master - release-ready version of the library
* develop - development version of the library
 
and could contain various supporting branches for new features and hotfixes. 

As a contributor you'll most likely be adding new features or fixing bugs in the development version of the library. 
This means that for each contribution you should create a new branch originating from the develop branch, 
modify it and send a pull request in order to merge it, again with the develop branch.

### Create new branch for your work

Make sure you're in develop branch running

    git branch -vv

you should see something like this: 

    * develop a76b4be [origin/develop] Update issue templates

Now you should pick **a name for your new branch that doesn't already exist**. 
The following returns a list of all the existing remote branches

    git branch -a

Alternatively, you can check them on `GitHub`.

Assume you want to add some new functionality. 
Then you have to create a new branch e.g. `feature/my_cool_new_feature`

Create new local branch

    git branch feature/my_cool_new_feature
    git checkout feature/my_cool_new_feature

push it to your fork

    git push -u my_fork feature/my_cool_new_feature

Note that the `-u` switch also sets up the tracking of the remote branch. 
Your new branch now is created!

### Verify your new branch (optional)

Now if you check the branches present on your repository
you'll see the `develop` and `master` branches as well as the one you just created

```bash
        user@mypc:~/dingo$git branch -vv
        develop        f82fcce [origin/develop] Revert "Revert "Update issue templates""
        * dingo_tutorial 1806b75 [origin/dingo_tutorial] notebook moved under /tutorials
        pairs          17d6d0b [origin/pairs] ignore notebook checkpoints
```

Note that without the `-u` switch you wouldn't see the tracking information for your new branch.

Your newly created remote branch is also available on GitHub
on your fork repository!

![branch_on_github](https://raw.githubusercontent.com/hariszaf/dingo/dingo_tutorial/tutorials/figs/branches_github.png)

Notice, we are **not** on the `dingo` repository under the `GeomScale` organization, but on the user's personal account. 

## Modify the branch (implement, implement, implement)

Before contributiong to a library by adding a new feature, or a bugfix, or improving documentation, 
it is always wise to interact with the community of developers, for example by opening an issue.

### Tests

Tests are placed in the `test` directory and use the [doctest](https://github.com/onqtam/doctest) library. 

It is recommended to add new test whenever you contribute a new functionality/feature.
Also if your contribution is a bugfix then consider adding this case to the test-suite.

### Push

At the end, push your changes to the remote branch

    git push my_fork feature/my_cool_new_feature

or if your local branch is tracking the remote one, just

    git push

## Pull request (the joy of sharing)

After pushing your work you should be able to see it on `GitHub`.

Click "Compare and pull request" button or the "New pull request" button.

Add title and description

![RP](https://raw.githubusercontent.com/hariszaf/dingo/dingo_tutorial/tutorials/figs/pr.png)

and click the "Create pull request" button.

## Review (ok this is not an exam)

After creating a pull request your code will be reviewed. You can propose one or more reviewers 
by clicking on the "Reviewers" button

![reviewer](https://user-images.githubusercontent.com/3660366/72349476-44ecc600-36e5-11ea-81cd-d0938d923529.png)

If there are no objections your changes will be merged. 
Otherwise you'll see some comments under the pull request and/or under specific lines of your code.
Then you have to make the required changes, commit them and push to your branch. 
Those changes will automatically be a part of the same pull request. This procedure will be repeated until the code 
is ready for merging.

If you're curious how it looks like you may see one of the open or closed 
[pull requests](https://github.com/GeomScale/dingo/pulls).