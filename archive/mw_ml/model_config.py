import numpy as np
from sklearn.ensemble import (
    RandomForestClassifier,
    GradientBoostingClassifier,
    AdaBoostClassifier,
    HistGradientBoostingClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from xgboost import XGBClassifier

MODEL_CONFIGS = {
    "logistic_regression": {
        "class": LogisticRegression,
        "fixed_params": {"max_iter": 1000},
        "search_space": {
            "C": [0.001, 0.01, 0.1, 1, 10, 100],
            "penalty": ["l2"],
            "solver": ["liblinear"],
        },
    },
    "random_forest": {
        "class": RandomForestClassifier,
        "fixed_params": {"random_state": 42},
        "search_space": {
            "n_estimators": [50, 100, 200],
            "max_depth": [None, 10, 20, 30],
            "min_samples_split": [2, 5, 10],
            "min_samples_leaf": [1, 2, 4],
            "bootstrap": [True, False],
        },
    },
    "svm": {
        "class": SVC,
        "fixed_params": {"probability": True},
        "search_space": {
            "C": [0.1, 1, 10, 100],
            "kernel": ["linear", "rbf", "poly"],
            "gamma": ["scale", "auto"],
        },
    },
    "knn": {
        "class": KNeighborsClassifier,
        "fixed_params": {},
        "search_space": {
            "n_neighbors": [3, 5, 7, 9],
            "weights": ["uniform", "distance"],
            "metric": ["euclidean", "manhattan"],
        },
    },
    "naive_bayes": {
        "class": GaussianNB,
        "fixed_params": {},
        "search_space": {"var_smoothing": np.logspace(-9, -1, 10)},
    },
    "xgboost": {
        "class": XGBClassifier,
        "fixed_params": {"random_state": 42, "eval_metric": "mlogloss"},
        "search_space": {
            "n_estimators": [50, 100, 200],
            "learning_rate": [0.0001, 0.001, 0.01, 0.1],
            "max_depth": [3, 5, 7, 10],
            "subsample": [0.6, 0.8, 1.0],
            "colsample_bytree": [0.6, 0.8, 1.0],
        },
    },
    "mlp": {
        "class": MLPClassifier,
        "fixed_params": {"max_iter": 1000, "random_state": 42},
        "search_space": {
            "hidden_layer_sizes": [(50,), (100,), (100, 50), (50, 100), (64, 128, 64)],
            "activation": ["tanh", "relu"],
            "solver": ["adam"],
            "alpha": [0.00001, 0.0001, 0.001, 0.01],
            "learning_rate": ["constant", "adaptive"],
            "learning_rate_init": [0.001, 0.0001, 0.00001],
        },
    },
    "gradient_boosting": {
        "class": GradientBoostingClassifier,
        "fixed_params": {"random_state": 42},
        "search_space": {
            "n_estimators": [50, 100, 200],
            "learning_rate": [0.1, 0.01, 0.001],
            "max_depth": [3, 5, 7],
            "min_samples_split": [2, 5, 10],
            "min_samples_leaf": [1, 2, 4],
        },
    },
    "hist_gradient_boosting": {
        "class": HistGradientBoostingClassifier,
        "fixed_params": {"random_state": 42},
        "search_space": {
            "learning_rate": [0.1, 0.01, 0.001],
            "l2_regularization": [0.0, 0.1, 0.01],
            "min_samples_leaf": [10, 15, 20],
        },
    },
    "adaboost": {
        "class": AdaBoostClassifier,
        "fixed_params": {"random_state": 42, "algorithm": "SAMME"},
        "search_space": {
            "n_estimators": [50, 100, 200],
            "learning_rate": [0.1, 0.01, 0.001],
        },
    },
    "qda": {
        "class": QuadraticDiscriminantAnalysis,
        "fixed_params": {},
        "search_space": {"reg_param": [0.0, 0.1], "tol": [1e-4, 1e-3]},
    },
}


def get_random_hyperparameters(model_name):
    """Helper function to generate random hyperparameters for a given model type"""
    search_space = MODEL_CONFIGS[model_name]["search_space"]
    random_params = {}
    for k, v in search_space.items():
        if isinstance(v[0], list) or isinstance(v[0], tuple):
            random_params[k] = v[np.random.randint(len(v))]
        else:
            random_params[k] = (
                str(np.random.choice(v))
                if isinstance(v[0], str)
                else np.random.choice(v)
            )
    return random_params
