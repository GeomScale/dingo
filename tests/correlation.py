
from dingo.illustrations import corr
from dingo import MetabolicNetwork

model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
corr(model)