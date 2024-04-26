class ReconstructionTimeNotSet(Exception):
    """raise this exception when the reconstruction time is None."""

    def __init__(self):
        super().__init__(
            "The reconstruction time has not been set yet. Set `PlotTopologies.time` before calling plotting functions."
        )
