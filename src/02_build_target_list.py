
import pandas as pd
from pathlib import Path

catalog = Path('data/catalog')

confirmed = pd.read_csv(catalog / 'confirmed_planets.csv')
confirmed = confirmed[confirmed['tic_id'].notna()].copy()
confirmed['label'] = 1

targets = confirmed[['tic_id','label']].drop_duplicates()
targets.to_csv(catalog / 'targets.csv', index=False)

print('targets.csv generado')
