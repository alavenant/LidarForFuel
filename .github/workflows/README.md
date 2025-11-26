# GitHub Actions Workflows

Ce dossier contient les workflows GitHub Actions pour l'intégration continue.

## Workflows disponibles

### test-simple.yml
Workflow simplifié qui exécute les tests sur Ubuntu avec Python 3.10.
- **Déclenchement** : Push et Pull Requests sur `main`, `master`, `develop`
- **OS** : Ubuntu latest
- **Python** : 3.10
- **Durée** : ~2-3 minutes

### test.yml
Workflow complet qui exécute les tests sur plusieurs OS et versions de Python.
- **Déclenchement** : Push et Pull Requests sur `main`, `master`, `develop`
- **OS** : Ubuntu latest, macOS latest
- **Python** : 3.8, 3.9, 3.10, 3.11
- **Fonctionnalités** : Tests avec couverture de code
- **Durée** : ~10-15 minutes

## Utilisation

Les workflows s'exécutent automatiquement à chaque push ou pull request.

Pour exécuter les tests localement :

```bash
# Activer l'environnement conda
conda activate lidarforfuel

# Exécuter les tests
pytest test/ -v

# Avec couverture
pytest test/ --cov=python --cov-report=html
```

## Configuration

Les workflows utilisent :
- `environment.yml` pour créer l'environnement conda
- `test/requirements.txt` pour les dépendances de test supplémentaires
- `pytest` pour exécuter les tests

