# Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>
#
# SPDX-License-Identifier: MIT


import itertools
import sys
from typing import Dict, List, Tuple, Union

import h5py
import hictkpy
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from sklearn.metrics import confusion_matrix

from . import IO


def average_within_threshold(numbers, threshold):
    i = 0
    while i < len(numbers) - 1:
        j = i + 1
        while j < len(numbers) and abs(numbers[j] - numbers[i]) <= threshold:
            j += 1
        if j > i + 1:
            average = sum(numbers[i:j]) / (j - i)
            numbers[i:j] = [average] * (j - i)
        i = j
    return numbers


def initialize_evaluation_dictionary() -> Dict[str, Union[List, int]]:
    """
    Initialize dictionary that contains relevant information to evaluate a caller.

    Returns
    -------
    Dict[str, Union[List, int]]
        dictionary mapping strings to either lists (0-1 vectors) or integers (counts)
    """
    keys1 = ["is_candidate_good", "is_anchor_found", "classification_vector", "GT_classification_vector"]
    keys2 = ["TN", "FP", "FN", "TP"]
    M = {**{key: [] for key in keys1}, **{key: 0 for key in keys2}}

    return M


def _is_anchor_in_stripes(GT_anchors: NDArray, pred_HIoIs: NDArray) -> Tuple[NDArray, NDArray]:
    """
    Check whether a ground truth anchor is inside a predicted stripe (i.e., in its horizontal domain).

    Parameters
    ----------
    GT_anchors: NDArray
        numpy array containing the coordinates of the ground truth anchors.
    pred_HIoIs: NDArray
        numpy array containing the horizontal domains

    Returns
    -------
    Tuple[NDArray, NDArray]
        two numpy arrays, expressing:

        * whether an anchor site is found (entry is set to 1) or not found (entry is set to 0)
        * whether a predicted stripe contains (entry is set to 1) or does not contain (entry is set to 0) an anchor
    """

    # Initialize arrays for is_anchor_found and is_candida_good with zeros
    is_anchor_found = np.zeros(len(GT_anchors), dtype=int)
    is_candida_good = np.zeros(len(pred_HIoIs), dtype=int)

    # Convert GT_anchors and pred_HIoIs to NumPy arrays for efficient comparison
    GT_anchors_arr = np.array(GT_anchors)
    pred_HIoIs_arr = np.array(pred_HIoIs)

    # Find if each GT anchor is contained in a candidate stripe
    for i, GT_anchor in enumerate(GT_anchors_arr):
        indices = (
            np.where(np.logical_and((GT_anchor >= pred_HIoIs_arr[:, 0]), (GT_anchor <= pred_HIoIs_arr[:, 1])))[0]
            if pred_HIoIs_arr.shape[0] > 0
            else []
        )

        if len(indices) > 0:
            is_anchor_found[i] = 1
            is_candida_good[indices] = 1

    return is_anchor_found, is_candida_good


def compute_classification_measures(
    confusion_matrices: Dict[str, Dict[str, int]], metrics: Dict[str, Dict[str, float]]
):
    """
    Compute classification measures out of confusion matrices, stored in dictionaries of dictionaries.

    Parameters
    ----------
    confusion_matrices: Dict[str, Dict[str, int]]
        dictionary of dictionaries, it stores a confusion matrix for each stripe caller.
    metrics: Dict[str, Dict[str, float]]
        dictionary of dictionaries, it stores classification measures for each stripe caller.
    """

    for prefix in confusion_matrices.keys():

        TP = confusion_matrices[prefix]["TP"]
        TN = confusion_matrices[prefix]["TN"]
        FP = confusion_matrices[prefix]["FP"]
        FN = confusion_matrices[prefix]["FN"]

        # Sensitivity, recall, hit rate, or true positive rate:
        TPR = TP / (TP + FN) if TP + FN > 0 else 0
        metrics[prefix]["TPR"] = TPR

        # Specificity, selectivity, or true negative rate:
        TNR = TN / (TN + FP) if TN + FP > 0 else 0
        metrics[prefix]["TNR"] = TNR

        # Precision or positive predictive value:
        PPV = TP / (TP + FP) if TP + FP > 0 else 0
        metrics[prefix]["PPV"] = PPV

        # Balanced accuracy:
        bACC = 0.5 * (TPR + TNR)
        metrics[prefix]["bACC"] = bACC

        # G-mean (i.e., geometric mean):
        GM = np.sqrt(TPR * TNR)
        metrics[prefix]["GM"] = GM

        # F1-score:
        F1c = 2 * (TPR * PPV) / (TPR + PPV) if TPR + PPV > 0 else 0
        metrics[prefix]["F1c"] = F1c

        # Fowlkes–Mallows index:
        FMc = np.sqrt(TPR * PPV)
        metrics[prefix]["FMc"] = FMc

        JI = TP / (TP + FN + FP)
        metrics[prefix]["JI"] = JI


def compute_recognition_measures(
    is_anchor_found: NDArray, is_candidate_good: NDArray, metrics: Dict[str, Dict[str, float]]
):
    """
    Compute classification measures out of confusion matrices, stored in dictionaries of dictionaries.

    Parameters
    ----------
    is_anchor_found: NDArray
        numpy array coding whether an anchor site is found (entry is set to 1) or not found (entry is set to 0)
    is_candidate_good: NDArray
        numpy array coding whether a predicted stripe contains (entry is set to 1) or does not contain (entry is set to
        0) an anchor
    metrics: Dict[str, Dict[str, float]]
        dictionary of dictionaries, it stores recognition measures for each stripe caller.
    """

    for prefix in metrics.keys():

        metrics[prefix]["AHR"] = np.sum(is_anchor_found[prefix]) / len(is_anchor_found[prefix])
        metrics[prefix]["FGC"] = (
            np.sum(is_candidate_good[prefix]) / len(is_candidate_good[prefix])
            if len(is_candidate_good[prefix]) > 0
            else 0
        )
        metrics[prefix]["F1r"] = (
            (2 * (metrics[prefix]["AHR"] * metrics[prefix]["FGC"]) / (metrics[prefix]["AHR"] + metrics[prefix]["FGC"]))
            if (metrics[prefix]["AHR"] + metrics[prefix]["FGC"]) > 0
            else 0
        )
        metrics[prefix]["FMr"] = np.sqrt(metrics[prefix]["AHR"] * metrics[prefix]["FGC"])


def _unique_HIoIs(HIoIs: NDArray) -> NDArray:
    """
    This function removes possible duplicates from the horizontal domains.

    Parameters
    ----------
    HIoIs: NDArray
       numpy array containing the predicted horizontal domains

    Returns
    -------
    NDArray
       numpy array containing the predicted horizontal domains without repetitions
    """

    uHIoIs = []
    seen = set()
    for HIoI in HIoIs:
        t_HIoI = tuple(HIoI)
        if t_HIoI not in seen:
            uHIoIs.append(HIoI)
            seen.add(t_HIoI)

    return np.array(uHIoIs)


def _compare_predictions_to_StripeBench_GT(
    GT: Dict[str, Dict[str, NDArray]], HIoIs: Dict[str, NDArray], clas_vec: Dict[str, NDArray]
) -> Dict[str, Union[List, int]]:
    """
    This function compares predictions of a given stripe caller against ground truth annotations.

    Parameters
    ----------
    GT: Dict[str, Dict[str, NDArray]]
       dictionary containing the ground truth annotations as defined in IO.py
    HIoIs: Dict[str, NDArray]
       dictionary containing the horizontal domains as defined in IO.py
    clas_vec: Dict[str, NDArray]
       dictionary containing the predicted classification vectors as defined in IO.py

    Returns
    -------
    Dict[str, Union[List, int]]
       dictionary mapping strings to either lists (0-1 vectors) or integers (counts)
    """

    # Remove repeating pairs (Chromosight and StripeCaller break stripes into sub-stripes):
    HIoIs["LT"] = _unique_HIoIs(HIoIs["LT"])
    HIoIs["UT"] = _unique_HIoIs(HIoIs["UT"])

    # Initialize dictionary that will contain information for the lower- upper-triangular matrices:
    M = dict()

    # Check which GT anchors fall into the candidate stripes AND which candidate stripes contain an anchor site:
    L_is_anchor_found, L_is_candidate_good = _is_anchor_in_stripes(
        GT_anchors=GT["anchors"]["LT"], pred_HIoIs=HIoIs["LT"]
    )
    U_is_anchor_found, U_is_candidate_good = _is_anchor_in_stripes(
        GT_anchors=GT["anchors"]["UT"], pred_HIoIs=HIoIs["UT"]
    )
    M["is_anchor_found"] = np.concatenate((L_is_anchor_found, U_is_anchor_found)).tolist()
    M["is_candidate_good"] = np.concatenate((L_is_candidate_good, U_is_candidate_good)).tolist()

    # Confusion matrix:
    y_true = GT["clas_vec"]["LT"].tolist() + GT["clas_vec"]["UT"].tolist()
    y_pred = np.concatenate((clas_vec["LT"], clas_vec["UT"]))
    M["TN"], M["FP"], M["FN"], M["TP"] = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    return M


def _compute_measures_StripeBench(M) -> Tuple[float]:
    """
    This function compares predictions of a given stripe caller against ground truth annotations.

    Parameters
    ----------
    M: Dict[str, Union[List, int]]
        dictionary mapping strings to either lists (0-1 vectors) or integers (counts)

    Returns
    -------
    Tuple[float]
        tuple that contains classification (TPR, TNR, PPV, bACC, GM, JI, F1c, FMc) and recognition (AHR, FGC, F1r, FMr)
        measures
    """

    # 1) STRIPE-BASED RECOGNITION MEASURES

    # Anchor Hit Rate:
    AHR = np.sum(M["is_anchor_found"]) / len(M["is_anchor_found"])

    # Fraction Good Candidates:
    FGC = np.sum(M["is_candidate_good"]) / len(M["is_candidate_good"]) if len(M["is_candidate_good"]) > 0 else 0

    # F1 and FM scores:
    F1r = 2 * (AHR * FGC) / (AHR + FGC) if AHR + FGC > 0 else 0
    F1r = F1r if ~np.isnan(F1r) else 0
    FMr = np.sqrt(AHR * FGC)

    # 2) ANCHOR-BASED CLASSIFICATION MEASURES

    # Sensitivity, recall, hit rate, or true positive rate:
    TPR = M["TP"] / (M["TP"] + M["FN"]) if M["TP"] + M["FN"] > 0 else 0

    # Specificity, selectivity, or true negative rate:
    TNR = M["TN"] / (M["TN"] + M["FP"]) if M["TN"] + M["FP"] > 0 else 0

    # Precision or positive predictive value:
    PPV = M["TP"] / (M["TP"] + M["FP"]) if M["TP"] + M["FP"] > 0 else 0

    # Balanced accuracy:
    bACC = 0.5 * (TPR + TNR)

    # G-mean (i.e., geometric mean):
    GM = np.sqrt(TPR * TNR)

    # Jaccard Index
    JI = M["TP"] / (M["TP"] + M["FP"] + M["FN"]) if (M["TP"] + M["FP"] + M["FN"]) > 0 else 0

    # F1 and FM scores:
    F1c = 2 * (TPR * PPV) / (TPR + PPV) if (TPR + PPV) > 0 else 0
    F1c = F1c if ~np.isnan(F1c) else 0
    FMc = np.sqrt(TPR * PPV)

    return TPR, TNR, PPV, bACC, GM, JI, AHR, FGC, F1c, FMc, F1r, FMr
