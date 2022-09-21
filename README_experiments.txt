
vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ git submodule update --init
vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ wget -O boost_1_76_0.tar.bz2 https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2
vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ tar xjf boost_1_76_0.tar.bz2
vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ rm boost_1_76_0.tar.bz2
vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ poetry shell
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ poetry install
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ pip3 install -i https://pypi.gurobi.com gurobipy
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ pip3 install PolyRound

(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ pip3 install hopsy

(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ mkdir simplified_transformed_polytopes
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ mkdir dingo_samples
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ mkdir dingo_samples_polyround_transformation
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ mkdir dingo_simplified_transformed_polytopes
(dingo-6kt0FIcQ-py3.8) vissarion@vissarion-ThinkPad-T490:~/workspace/dingo_clean$ python tests/experiments.py iEC1344_C 1
Set parameter Username
Academic license - for non-commercial use only - expires 2023-07-08
Simplify
Sampling
phase 1: number of correlated samples = 2100, effective sample size = 13
 [1]maximum marginal PSRF: 1.20919
Segmentation fault (core dumped)


