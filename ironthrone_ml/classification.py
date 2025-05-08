"""v1.2 - 241226"""

import os
import pandas as pd
import numpy as np
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    f1_score,
    precision_score,
    accuracy_score,
    make_scorer,
)
import numpy as np
import joblib
import re
from mw_ml import logging
from mw_ml.model_config import MODEL_CONFIGS, get_random_hyperparameters
from sklearn.ensemble import VotingClassifier
from collections import Counter


def custom_performance_metric(y_true, y_pred, y_value_of_interest, alpha=0.9):
    """
    Custom performance metric that combines overall F1 score and recall - false positive rate (FPR) for a specific class.

    Parameters:
    - y_true (array-like): True labels.
    - y_pred (array-like): Predicted labels.
    - y_value_of_interest (int or str): The class of interest.

    Returns:
    - float: The composite performance metric value.
    """
    # # * Combined Overall F1 score and Recall&FPR for a specific class
    # # Overall F1 score
    # overall_f1 = f1_score(y_true, y_pred, average="micro")

    # # Custom metric for y_value_of_interest
    # y_true_binary = (y_true == y_value_of_interest).astype(int)
    # y_pred_binary = (y_pred == y_value_of_interest).astype(int)
    # tn, fp, fn, tp = confusion_matrix(y_true_binary, y_pred_binary).ravel()
    # fpr = fp / (fp + tn)
    # recall = tp / (tp + fn)
    # custom_metric = recall - fpr

    # # Composite score
    # performance = overall_f1 + custom_metric

    # * Weighted combination of precision and recall (more emphasis on precision) of the y_value_of_interest
    report = classification_report(y_true, y_pred, output_dict=True, zero_division=0)
    precision = report[str(y_value_of_interest)]["precision"]
    recall = report[str(y_value_of_interest)]["recall"]
    performance = alpha * precision + (1 - alpha) * recall

    return performance


def train_classification_model(
    train_df,
    test_df,
    resultdir,
    target_column="Y",
    y_value_of_interest=None,
    models_to_run=[
        "logistic_regression",
        "random_forest",
        "knn",
        "naive_bayes",
        "xgboost",
        "mlp",
        "gradient_boosting",
        "hist_gradient_boosting",
        "adaboost",
        "qda",
    ],
    ensemble_hard=False,
    ensemble_soft=False,
    random_seed=1024,
):
    """
    Trains multiple classification models and selects the best one based on a custom performance metric or F1 score.
    Optionally trains ensemble models using hard and/or soft voting.

    Parameters:
    - train_df (DataFrame): Training dataset.
    - test_df (DataFrame): Testing dataset.
    - resultdir (str): Directory to save results and logs.
    - target_column (str): The name of the target column.
    - random_seed (int): Random seed for reproducibility.
    - y_value_of_interest (int or str, optional): The class of interest for the custom performance metric.
    - models_to_run (list, optional): List of model names to run. Defaults to all models.
    - ensemble_hard (bool, optional): Whether to train an ensemble model with hard voting. Defaults to False.
    - ensemble_soft (bool, optional): Whether to train an ensemble model with soft voting. Defaults to False.

    Returns:
    - tuple: y_test, y_pred_best, y_pred_prob_best, best_clf, best_models, best_model_name
    """

    # Separate features and target from both datasets
    X_train = train_df.drop(columns=[target_column])
    y_train = train_df[target_column]
    X_test = test_df.drop(columns=[target_column])
    y_test = test_df[target_column]

    # Data preprocessing (scaling)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Step 1: Define models and hyperparameter grids
    models_params = {
        model_name: {
            "model": MODEL_CONFIGS[model_name]["class"](
                **MODEL_CONFIGS[model_name]["fixed_params"]
            ),
            "params": MODEL_CONFIGS[model_name]["search_space"],
        }
        for model_name in models_to_run
    }

    # Filter models to run based on the provided list
    if models_to_run:
        models_params = {k: v for k, v in models_params.items() if k in models_to_run}

    # Define custom scorer
    if y_value_of_interest is not None:
        custom_scorer = make_scorer(
            custom_performance_metric,
            y_value_of_interest=y_value_of_interest,
            greater_is_better=True,
        )
    else:
        custom_scorer = "f1_micro"

    # Step 2: Loop through models, perform RandomizedSearchCV, and evaluate
    fitted_models = {}
    best_performance = -999
    best_clf = None
    y_pred_best = None
    best_model_name = None

    for model_name, mp in models_params.items():
        print(f"\nRunning hyperparameter search for {model_name}...")
        model = mp["model"]
        params = mp["params"]

        # Initialize RandomizedSearchCV
        search = RandomizedSearchCV(
            model,
            params,
            cv=5,
            n_iter=20,
            scoring=custom_scorer,
            refit=True,
            n_jobs=-1,
            random_state=random_seed,
            verbose=0,
        )
        search.fit(X_train_scaled, y_train)

        # Get the best model
        best_model = search.best_estimator_
        fitted_models[model_name] = best_model

        # Predict on test data and evaluate
        y_pred = best_model.predict(X_test_scaled)
        y_pred_prob = best_model.predict_proba(X_test_scaled)

        if y_value_of_interest is not None:
            performance = custom_performance_metric(y_test, y_pred, y_value_of_interest)
        else:
            performance = f1_score(y_test, y_pred, average="micro")

        logging(
            f"{model_name} - : Performance: {performance:.4f}",
            resultdir,
            f"log.txt",
        )
        logging(
            f"Classification Report for {model_name}:\n{classification_report(y_test, y_pred)}",
            resultdir,
            f"log.txt",
        )

        if performance > best_performance:
            best_performance = performance
            best_clf = best_model
            y_pred_best = y_pred
            y_pred_prob_best = y_pred_prob
            best_model_name = model_name

        # Save prediction results to the test dataframe
        test_df[f"{model_name}_pred"] = y_pred
        for i, class_name in enumerate(best_model.classes_):
            test_df[f"{model_name}_prob_{class_name}"] = y_pred_prob[:, i]

    test_df[f"{best_model_name}_pred_best"] = y_pred_best
    for i, class_name in enumerate(best_clf.classes_):
        test_df[f"{best_model_name}_prob_best_{class_name}"] = y_pred_prob_best[:, i]

    # Step 3: Output the best model and its predictions
    logging("\nBest of the Best Model: %s" % best_clf, resultdir, f"log.txt")
    logging("Best Model Test: %s" % best_performance, resultdir, f"log.txt")
    logging(
        f"Classification Report for Best Model {best_model_name}:\n{classification_report(y_test, y_pred_best)}",
        resultdir,
        f"log.txt",
    )

    # === Save the best model === #
    best_params = search.best_params_
    # Create a string representation of the hyperparameters
    params_str = "_".join(f"{key}={value}" for key, value in best_params.items())
    # Replace any characters that are not safe for filenames
    safe_params_str = re.sub(r"[^A-Za-z0-9_=-]", "", params_str)
    filename = f"best_model_{best_model_name}_{safe_params_str}.pkl"
    joblib.dump(best_clf, os.path.join(resultdir, filename))

    # Train ensemble models if specified
    if ensemble_hard:
        print("\nTraining ensemble model with hard voting...")
        (
            y_test,
            y_pred_hard,
            _y_pred_prob_hard,
            confidence_scores_hard,
            _voting_clf_hard,
        ) = train_ensemble_voting_classifier(
            X_train_scaled,
            y_train,
            X_test_scaled,
            y_test,
            resultdir,
            use_soft_voting=False,
            models_to_run=models_to_run,
        )
        test_df["ensemble_hard_pred"] = y_pred_hard
        test_df["ensemble_hard_confidence"] = confidence_scores_hard
        fitted_models["ensemble_hard"] = _voting_clf_hard
        # Save the ensemble model
        joblib.dump(
            _voting_clf_hard, os.path.join(resultdir, "ensemble_model_hard.pkl")
        )

    if ensemble_soft:
        print("\nTraining ensemble model with soft voting...")
        (
            y_test,
            y_pred_soft,
            y_pred_prob_soft,
            _confidence_scores_soft,
            voting_clf_soft,
        ) = train_ensemble_voting_classifier(
            X_train_scaled,
            y_train,
            X_test_scaled,
            y_test,
            resultdir,
            use_soft_voting=True,
            models_to_run=models_to_run,
        )
        test_df["ensemble_soft_pred"] = y_pred_soft
        for i, class_name in enumerate(voting_clf_soft.classes_):
            test_df[f"ensemble_soft_prob_{class_name}"] = y_pred_prob_soft[:, i]
        fitted_models["ensemble_soft"] = voting_clf_soft
        # Save the ensemble model
        joblib.dump(voting_clf_soft, os.path.join(resultdir, "ensemble_model_soft.pkl"))

    return (
        test_df,
        fitted_models,
        best_model_name,
        y_pred_best,
        y_pred_prob_best,
        best_clf,
    )


# ============================ #
# ====== ENSEMBLE MODEL ====== #
# ============================ #
from sklearn.ensemble import VotingClassifier
from collections import Counter


def train_ensemble_voting_classifier(
    X_train_scaled,
    y_train,
    X_test_scaled,
    y_test,
    resultdir,
    use_soft_voting=False,
    models_to_run=[
        "logistic_regression",
        "random_forest",
        "knn",
        "naive_bayes",
        "xgboost",
        "mlp",
        "gradient_boosting",
        "hist_gradient_boosting",
        "adaboost",
    ],
    n_models=3,
):
    """
    Trains an ensemble voting classifier using multiple models with random hyperparameters.

    Parameters:
    - X_train_scaled (array-like): Scaled training features.
    - y_train (array-like): Training labels.
    - X_test_scaled (array-like): Scaled testing features.
    - y_test (array-like): Testing labels.
    - resultdir (str): Directory to save results and logs.
    - use_soft_voting (bool): Whether to use soft voting. Defaults to False.
    - random_seed (int): Random seed for reproducibility.
    - models_to_run (list): List of model types to include in the ensemble
    - n_models (int): Number of variations to create for each model type

    Returns:
    - tuple: y_test, y_pred, y_pred_prob, confidence_scores, voting_clf
    """
    # Initialize empty model list
    model_list = []

    # Generate models for each selected type
    for model_name in models_to_run:
        if model_name not in MODEL_CONFIGS:
            print(f"Warning: {model_name} is not a supported model type")
            continue

        model_class = MODEL_CONFIGS[model_name]["class"]

        for i in range(n_models):
            # Get random hyperparameters
            params = get_random_hyperparameters(model_name)
            params.update(MODEL_CONFIGS[model_name]["fixed_params"])

            # Create model instance
            model = model_class(**params)
            model_list.append((f"{model_name}_{i + 1}", model))

    # Initialize VotingClassifier
    voting_type = "soft" if use_soft_voting else "hard"
    voting_clf = VotingClassifier(estimators=model_list, voting=voting_type, n_jobs=-1)

    # Train the VotingClassifier on the scaled training data
    voting_clf.fit(X_train_scaled, y_train)

    # Make predictions
    y_pred = voting_clf.predict(X_test_scaled)

    # If using soft voting, also get predicted probabilities
    y_pred_prob = None
    if use_soft_voting:
        y_pred_prob = voting_clf.predict_proba(X_test_scaled)

    # Calculate confidence scores for both voting types
    confidence_scores = []
    if use_soft_voting:
        # Confidence score same as the probability for the predicted class (Maximum among probabilities assigned to each class)
        for i in range(len(X_test_scaled)):
            pred_class = y_pred[i]
            confidence = y_pred_prob[i][
                pred_class
            ]  # Probability of the predicted class
            confidence_scores.append(confidence)
    else:
        # For hard voting, calculate how many models agreed on the prediction
        for i in range(len(X_test_scaled)):
            preds = [
                model.predict(X_test_scaled[i].reshape(1, -1)).item()
                for model in voting_clf.estimators_
            ]
            vote_counts = Counter(preds)
            ensemble_pred = vote_counts.most_common(1)[0][0]
            confidence = vote_counts[ensemble_pred] / len(voting_clf.estimators_)
            confidence_scores.append(confidence)

    # Evaluate the ensemble model
    logging(
        f"Voting Ensemble Accuracy: {accuracy_score(y_test, y_pred):.4f}",
        resultdir,
        f"log.txt",
    )
    logging(
        f"Voting Ensemble Classification Report:\n{classification_report(y_test, y_pred, zero_division=1)}",
        resultdir,
        f"log.txt",
    )

    return y_test, y_pred, y_pred_prob, confidence_scores, voting_clf
