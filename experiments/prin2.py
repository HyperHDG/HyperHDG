#!/usr/bin/env python3

import logging, os, subprocess, datetime, sys

class Logger(logging.Logger):
  """prin2.Logger
  use via:
     logging.setLoggerClass(prin2.Logger)
  """

  def __init__(self, name=os.path.basename(__file__), level=logging.INFO, log_dir="logs", file_mode="w"):
    super().__init__(name, level)

    self.log_dir = log_dir

    self.time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    self.git_hash =  subprocess.run(["git", "-C", os.path.dirname(__file__), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
    self.git_mod  = subprocess.run(["git", "-C", os.path.dirname(__file__), "status", "--porcelain", __file__], check=True, capture_output=True, text=True).stdout.strip()
    if self.git_mod:
      self.git_hash += "-dirty"

    log_formatter  = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

    os.system("mkdir -p logs")
    fhandler = logging.FileHandler(f"{self.log_dir}/{name}.{self.time_stamp}.log", mode=file_mode)
    fhandler.setFormatter(log_formatter)
    fhandler.setLevel(level)

    shandler = logging.StreamHandler(sys.stdout)
    shandler.setFormatter(log_formatter)
    shandler.setLevel(level)

    self.addHandler(fhandler)
    self.addHandler(shandler)

  def log_args(self, args):
    self.info(f"{self.name} started with")
    self.info(f"  git_hash={self.git_hash}")
    for key, val in vars(args).items():
      self.info(f"  {key}={val}")
