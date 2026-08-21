"""Public imports for sequential relative-volume updating.

Existing scripts imported ``dingo.volume_updating``.  Keep that import path,
but expose only the sequential wrapper around the supplied static routine.
"""

from dingo.dynamic_volume import (
    StaticVolumeResult,
    StaticVolumeUpdater,
    HalfspaceCut,
    MatlabVolumeUpdater,
    cuts_from_bounds,
)

__all__ = [
    "StaticVolumeResult",
    "StaticVolumeUpdater",
    "HalfspaceCut",
    "MatlabVolumeUpdater",
    "cuts_from_bounds",
]
