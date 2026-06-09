'''Exceptions raised by the superelliptic package.'''


class PendingMigrationError(NotImplementedError):
    '''Raised when a method depends on a legacy component that is not migrated yet.'''
