class NoHomologousSequenceException(Exception):
    pass


class AmbiguousQuerySequenceException(Exception):
    pass


class SubseqMatchError(Exception):
    pass


class UnknownStrand(Exception):
    pass


class IncorrectDatabaseChoice(Exception):
    def __init__(self):
        super.__init__("Incorrect database choice. This should never be reached.")


class SubprocessException(Exception):
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors


class RNAfoldException(SubprocessException):
    pass


class TurboFoldException(SubprocessException):
    pass


class LocarnaException(SubprocessException):
    pass


class MuscleException(SubprocessException):
    pass


class CentroidHomfoldException(SubprocessException):
    pass


class ClustaloException(SubprocessException):
    pass


class RNAalifoldException(SubprocessException):
    pass


class RefoldException(SubprocessException):
    pass


class RNAfold(SubprocessException):
    pass


class HybridssminException(SubprocessException):
    pass


class CmemitException(SubprocessException):
    pass


class CmscanException(SubprocessException):
    pass


class CmalignException(SubprocessException):
    pass