'''Exceptions raised by the superelliptic package.'''


class PendingMigrationError(NotImplementedError):
    '''Raised when an expected migrated component cannot be loaded.'''
