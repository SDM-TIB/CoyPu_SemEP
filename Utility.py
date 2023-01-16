import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from rdflib import Graph
import scipy

from pykeen.triples import TriplesFactory
import pykeen
from pykeen.pipeline import pipeline

def convert_to_category(world_bank, indicator, lower):
    labels = ['low', 'mediumLow', 'medium', 'mediumHigh', 'high']
    if lower:
        labels.reverse()
    df = world_bank.loc[world_bank.indicator==indicator]
    df["value"] = pd.to_numeric(df["value"], downcast="float")
    df = df.sort_values(by='value')  # , ascending=boolean
    a = list(df['value'].values)
    n_split = np.array_split(a, 5)
    category = pd.cut(df.value, bins=[min(n_split[0])-1, max(n_split[0]), max(n_split[1]), max(n_split[2]),
                                      max(n_split[3]), max(n_split[4])], labels=labels)
    df.insert(3, 'category_indicator',category)
    df.drop(columns='value', inplace=True)
    return df


def load_dataset(path, name):
    triple_data = open(path + name).read().strip()
    data = np.array([triple.split(',') for triple in triple_data.split('\n')])
    tf_data = TriplesFactory.from_labeled_triples(triples=data)
    return tf_data, triple_data


def create_model(tf_training, tf_testing, embedding, n_epoch, path):
    results = pipeline(
        training=tf_training,
        testing=tf_testing,
        model=embedding,
        training_loop='sLCWA',
        #         negative_sampler='bernoulli',
        negative_sampler_kwargs=dict(
            filtered=True,
        ),
        # Training configuration
        training_kwargs=dict(
            num_epochs=n_epoch,
            use_tqdm_batch=False,
        ),
        # Runtime configuration
        random_seed=1235,
        device='gpu',
    )
    model = results.model
    results.save_to_directory(path + embedding)
    return model, results


def get_learned_embeddings(model):
    entity_representation_modules: List['pykeen.nn.RepresentationModule'] = model.entity_representations
    relation_representation_modules: List['pykeen.nn.RepresentationModule'] = model.relation_representations

    entity_embeddings: pykeen.nn.Embedding = entity_representation_modules[0]
    relation_embeddings: pykeen.nn.Embedding = relation_representation_modules[0]

    entity_embedding_tensor: torch.FloatTensor = entity_embeddings()
    relation_embedding_tensor: torch.FloatTensor = relation_embeddings()
    return entity_embedding_tensor, relation_embedding_tensor


def create_dataframe_predicted_entities(entity_embedding_tensor, entity, training):
    df = pd.DataFrame(entity_embedding_tensor.cpu().detach().numpy())
    df['target'] = list(training.entity_to_id)
    new_df = df.loc[df.target.isin(list(entity))]
    return new_df.iloc[:, :-1], new_df, df


def cosine_sim(x, y):
    return abs(1 - scipy.spatial.distance.cosine(x, y))


def matrix_similarity(new_df, f_dist, th):
    array = new_df.set_index('target')
    entity = list(array.index.values)
    sim_matrix = pd.DataFrame(index=entity, columns=entity)
    sim_matrix = sim_matrix.fillna(0.0)
    list_sim = []

    for index, row in array.iterrows():
        for indexC, rowC in array.iterrows():
            sim = f_dist(row.values, rowC.values)
            sim = round(sim, 5)
            sim_matrix.at[index, indexC] = sim
            list_sim.append(sim)

    threshold = np.percentile(list_sim, th)
    print("percentil", threshold)
    # for col in sim_matrix.columns:
    #    sim_matrix.loc[sim_matrix[col] < threshold, [col]] = 0
    return sim_matrix, threshold


# === Save cosine similarity matrix with the structure SemEP need
def SemEP_structure(name, sim_matrix, sep):
    f = open(name, mode="w+")
    f.write(str(sim_matrix.shape[0]) + "\n")
    f.close()
    sim_matrix.to_csv(name, mode='a', sep=sep, index=False, header=False, float_format='%.5f')


def create_entitie(list_n, ENTITIES_FILE):
    entities = "\n".join(str(x) for x in list_n)
    n_ent = str(len(list_n))
    entity = open(ENTITIES_FILE, mode="w+")
    entity.write(n_ent + "\n" + entities)
    entity.close()