# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [1.1.0] - 2026-05-28
### Added
- Created new `core` directory to handle shared logic.

### Changed
- **Major Architectural Refactor**: Transitioned codebase to a clean Model-View-Controller (MVC) structure.
- Extracted all FastF1 data loading and Pandas telemetry processing out of the UI layer into `core/data_loader.py`.
- Consolidated shared formatting and caching helpers into `core/utils.py`.
- Significantly reduced line counts in `race_replayer.py` and `analytics_dashboard.py`, separating frontend visualization from backend data crunching.

## [1.0.0] - 2026-04-15
### Added
- Initial stable release of F1 Visualizer Hub.
- Main launcher for both Race Replayer and Analytics Dashboard.
- Sidecar VisPy 3D rendering mode with interactive camera and HUD.
- Race playback controls, leaderboard, tyre/status indicators, and telemetry overlays.
- Offline cache support for FastF1 session data.

### Changed
- Refined sidecar 3D driver label rendering for improved readability and cleaner visuals.
