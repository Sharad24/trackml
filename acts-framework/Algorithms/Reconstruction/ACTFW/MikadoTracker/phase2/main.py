#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import time

import numpy as np

from trackml.dataset import load_dataset
from trackml.score import score_event

from model import Model


PATH_TO_DATA = '/home/data'
MAX_TIME_PER_EVENT = 600
ACCURACY_INF = 0.5
get_clock = time.perf_counter


def mixt_score(accuracy_mean, time_per_event):
    if time_per_event <= 0 or time_per_event > MAX_TIME_PER_EVENT:
        return -1  # Something went terribly wrong
    if accuracy_mean < ACCURACY_INF:
        return 0
    speed = np.log(1.0 + ( MAX_TIME_PER_EVENT / time_per_event ))
    score = np.sqrt(speed * (accuracy_mean - ACCURACY_INF)**2)
    return score


def main():
    tracker = Model()
    time_spent = 0
    n_event = 0
    score_sum = 0
    for event_id, hits, cells, truth in load_dataset(PATH_TO_DATA, parts=['hits', 'cells', 'truth']):
        print("Runing event", event_id, "...", flush=True)
        # Make predictions
        t_start = get_clock()
        sub = tracker.predict_one_event(event_id, hits, cells=cells)
        t_end = get_clock()
        # Compute accuracy score
        score = score_event(truth, sub)
        # accumulate time, score, number of events
        time_spent += t_end - t_start
        score_sum  += score
        n_event += 1
        time_per_event = time_spent / n_event
        score = score_sum / n_event
        # Print information
        print("event", event_id, "accuracy score :", score)
        print("event", event_id, 'time spent     :', t_end - t_start)
        print('total time spent:', time_spent)
        print("running speed   : {:.3f} sec per event".format(time_spent  / n_event))
        print("running score   :", mixt_score(score, time_per_event))
        print('-----------------------------------', flush=True)
        if n_event>100:
            break
    if n_event == 0:
        print("Warning: no event where found in the given directory.")
        exit()
    if time_spent <= 0:
        print("Warning : execution time <= 0. Something went wrong !")

    time_per_event = time_spent / n_event
    score = score_sum / n_event

    print("Accuracy mean      :", score)
    print("Time per event     :", time_per_event)
    print("Overall mixt score :", mixt_score(score, time_per_event))


if __name__ == '__main__':
    main()
