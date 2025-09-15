# SPDX-FileCopyrightText: Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

import time
import logging
logging.basicConfig()
from torch.multiprocessing import get_context

class BaseTask:
    def __init__(self, name):
        self.name = name

    def _process_start(self, ctx, params, in_q, out_q, exit_event):
        self.process = ctx.Process(target=self._process_loop, args=(params, in_q, out_q, exit_event))
        self.process.start()

    def _process_join(self):
        self.process.join()

    def _process_loop(self, params, in_q, out_q, exit_event):
        logging.info(f"{self.name}: process_func started")
        self.setup(params)
        while True:
            item = in_q.get()
            if not item: # None-like signals end of processing; forward it through the chain 
                out_q.put(item)
                break
            else:
                out_q.put(self._process_one(item))
        exit_event.wait()
        self.teardown()
        logging.info(f"{self.name}: process_func completed")

    def _process_one(self, item):
        logging.info(f"{self.name}: {item['tag']} started")
        start_time = time.time()
        try:
            item = self.process_one(item)
        except Exception as e:
            logging.error(f"{self.name}: {item['tag']}: ")
            logging.error(e)
        end_time = time.time()
        logging.info(f"{self.name}: {item['tag']} completed in {(end_time - start_time):0.2f} seconds")
        return item

    # Should be implemented by derived task
    def setup(self, params):
        return

    # Should be implemented by derived task
    def teardown(self):
        return

    # Should be implemented by derived task
    def process_one(self, item):
        return item

class BasePipeline:
    def __new__(cls, multiprocessing=False, tasks=[]):
        if cls is BasePipeline:
            if multiprocessing:
                return BasePipelineMultiprocessing(multiprocessing, tasks)
            else:
                return BasePipelineSequential(multiprocessing, tasks)
        else:
            return super().__new__(cls)

class BasePipelineMultiprocessing(BasePipeline):
    def __init__(self, multiprocessing=True, tasks=[]):
        self.ctx = get_context("spawn")
        self.tasks = tasks # a list of BaseTask instances, to be connected in series
        self.queues = [] # tasks[i] reads from queues[i] and writes to queues[i+1]
        self.exit_events = [] # Each task should terminate only after the corresponding exit_event is signalled
        for i in range(len(self.tasks)+1):
            self.queues.append(self.ctx.Queue())
        for i in range(len(self.tasks)):
            self.exit_events.append(self.ctx.Event())

    def start(self, params):
        logging.info(f"pipeline: Starting in multiprocessing mode")
        for i in range(len(self.tasks)):
            self.tasks[i]._process_start(self.ctx, params, self.queues[i], self.queues[i+1], self.exit_events[i])

    def run(self, inputs):
        for item in inputs:
            logging.info(f"pipeline: Sending {item['tag']}")
            self.queues[0].put(item)
        outputs = []
        for item in inputs:
            outputs.append(self.queues[-1].get())
            logging.info(f"pipeline: Received {item['tag']}")
        return outputs

    def end(self):
        for event in self.exit_events:
            event.set()
        self.queues[0].put(None) # Signals exit to all tasks 
        for task in self.tasks:
            task._process_join()
        for queue in self.queues:
            queue.close()
        logging.info(f"pipeline: Completing run")

class BasePipelineSequential(BasePipeline):
    def __init__(self, multiprocessing=False, tasks=[]):
        self.tasks = tasks

    def start(self, params):
        logging.info(f"pipeline: Starting in sequential mode")
        for task in self.tasks:
            task.setup(params)

    def run(self, inputs):
        for item in inputs:
            logging.info(f"pipeline: Sending {item['tag']}")
            for task in self.tasks:
                item = task._process_one(item)
            logging.info(f"pipeline: Received {item['tag']}")

    def end(self):
        for task in self.tasks:
            task.teardown()
        logging.info(f"pipeline: Completing run")
