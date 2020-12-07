import os
import pysam
import gzip
import sqlite3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_genome(root_dir, accession, start, end):
    ff = resolver(root_dir, accession) + '.fasta.gz'

    g = pysam.FastaFile(ff)
    extracted = g.fetch(reference=accession, start=start, end=end)
    g.close()

    ace = accession.encode()
    with gzip.open(ff) as gz:
        for line in gz:
            if line[1:].startswith(ace):
                return SeqRecord(Seq(extracted), id=accession, description=line[1:].decode().strip())

    raise RuntimeError("Failed to find accession {} in {}".format(accession, ff))


def resolver(root, accession):
    """Return path for accession

    Path for genome file. Includes subdirectories based on last 4 characters from accession.
    The subdirectories are created if not exists

    :param root: PATH where the genome database should be created
    :param accession: genome pointer based on ACCESSION.VERSION
    :return:
    """

    ac = accession.split('.')[0]
    acc = ac.replace('_', '')
    dp = os.path.join(root, acc[:2], acc[2:4], acc)
    return os.path.join(dp, accession)


class GenomeDB(object):
    def __init__(self, dbfile, dbroot=None):
        self.con = sqlite3.connect(dbfile, timeout=300)
        self.con.row_factory = sqlite3.Row

        if dbroot:
            self.dbroot = dbroot
        else:
            self.dbroot = os.path.dirname(dbfile)

        self.tmpdir = os.path.join(self.dbroot, '__tmpdir')
        os.makedirs(self.tmpdir, exist_ok=True)

    def _read_db(self, accession):
        with self.con:
            c = self.con.execute("SELECT description, seq, gff FROM genomes WHERE acc = ?;", (accession,))
            return c.fetchone()

    def load(self, accession, start=None, end=None):

        row = self._read_db(accession)

        if row is None:
            raise KeyError("Accession {} is not avalible in the database.".format(accession))
        else:
            return LoadSeq(
                dbroot=self.dbroot,
                accession=accession,
                start=start,
                end=end,
                dbseq=self._decompress_or_none(row["seq"]),
                dbgff=self._decompress_or_none(row["gff"]),
                dbdescription=row["description"]
            )

    @staticmethod
    def _decompress_or_none(data):
        if data is None:
            return None
        else:
            return gzip.decompress(data).decode()

    def load_genome(self, accession, start, end):
        data = self.load(accession, start, end)
        return data.load_seq()

    def load_gff(self, accession):
        data = self.load(accession)
        return data.load_gff()

    def load_genome_and_gff(self, accession, start, end):
        data = self.load(accession, start, end)
        return data.load_seq(), data.load_gff()

    def load_db(self):
        with self.con:
            c = self.con.execute("SELECT acc FROM genomes;")
            return {i[0] for i in c}

    def acc_in_db(self, accession: str) -> bool:
        with self.con:
            c = self.con.execute("SELECT acc FROM genomes WHERE acc = ?;", (accession,))
            if c.fetchone() is None:
                return False
            else:
                return True

    @staticmethod
    def _tuplify(x):
        if not isinstance(x, tuple):
            return (x,)
        else:
            return x

    def missing_accessions(self, accessions: (set, list)) -> set:
        with self.con:
            try:
                self.con.execute("CREATE TEMPORARY TABLE temp(acc TEXT PRIMARY KEY NOT NULL);")
                self.con.executemany("INSERT INTO temp VALUES(?);", map(self._tuplify, accessions))
                c = self.con.execute(
                    "SELECT temp.acc FROM temp LEFT JOIN genomes ON temp.acc = genomes.acc WHERE genomes.acc IS NULL;"
                )

                return {r[0] for r in c}

            finally:
                self.con.execute("DROP TABLE temp;")

    def accessions_in_db(self, accessions: (set, list)) -> set:
        with self.con:
            try:
                self.con.execute("CREATE TEMPORARY TABLE temp(acc TEXT PRIMARY KEY NOT NULL);")
                self.con.executemany("INSERT INTO temp VALUES(?);", map(self._tuplify, accessions))
                c = self.con.execute(
                    "SELECT temp.acc FROM temp INNER JOIN genomes ON temp.acc = genomes.acc;"
                )

                return {r[0] for r in c}

            finally:
                self.con.execute("DROP TABLE temp;")

    def close(self):
        self.con.close()


class LoadSeq(object):
    def __init__(self, dbroot, accession, start, end, dbseq, dbgff, dbdescription):
        self.dbroot = dbroot
        self.accession = accession
        self.start = start
        self.end = end
        self.dbseq = dbseq
        self.dbgff = dbgff
        self.dbdescription = dbdescription

    def load_seq(self):
        if self.dbseq is None:
            return load_genome(self.dbroot, self.accession, self.start, self.end)

        else:
            if self.start is None:
                _seq = self.dbseq
            else:
                _seq = self.dbseq[self.start:self.end]

            return SeqRecord(
                Seq(_seq),
                id=self.accession,
                description=self.dbdescription
            )

    def load_gff(self):
        # do if on "dbseq" rather then dbgff
        # the dbgff might be empty even if seq is in db (no annotation on sequence - common for short sequences)

        if self.dbseq is None:
            # load gff from hdd
            raise NotImplementedError
        else:
            # add gff headers and return string
            if self.dbgff is None:
                return None
            else:
                gffstring = "##gff-version 3\n##sequence-region {} 1 {}\n".format(self.accession, len(self.dbseq)) \
                    + self.dbgff + "\n###"

                return gffstring