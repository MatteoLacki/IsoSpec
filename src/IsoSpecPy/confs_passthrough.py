



class ConfsPassthrough(object):
    def __init__(self, confs_parser, size):
        self.confs_parser = confs_parser
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        return self.confs_parser(idx)


