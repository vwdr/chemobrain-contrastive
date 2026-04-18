.PHONY: help install test smoke lint clean download preprocess train analyze

help:
	@echo "Targets:"
	@echo "  install     - pip install -r requirements.txt"
	@echo "  smoke       - run the synthetic-data smoke test"
	@echo "  test        - run pytest on tests/"
	@echo "  download    - fetch confirmed GEO datasets"
	@echo "  preprocess  - QC + merge + HVG + cell-type annotation"
	@echo "  train       - train MC-ContrastiveVI (needs merged data)"
	@echo "  analyze     - gene attribution + variance decomposition on latest run"
	@echo "  clean       - remove __pycache__, *.pyc"

install:
	pip install -r requirements.txt

smoke:
	python tests/test_smoke.py

test:
	pytest tests/ -v

download:
	python scripts/00_download_data.py

preprocess:
	python scripts/01_preprocess.py

train:
	python scripts/02_train.py --config configs/default.yaml

analyze:
	@latest=$$(ls -1dt runs/*/ | head -n1); \
	echo "Analyzing $$latest"; \
	python scripts/03_analyze.py --run $$latest

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete
