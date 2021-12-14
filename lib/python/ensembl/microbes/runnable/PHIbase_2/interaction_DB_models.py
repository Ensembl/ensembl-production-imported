# coding: utf-8
from sqlalchemy import (
        Column,
        DECIMAL,
        DateTime,
        TIMESTAMP,
        Enum,
        Float,
        ForeignKey,
        Index,
        Identity,
        String,
        Table,
        Text,
        text,
        )
from sqlalchemy.dialects.mysql import (
        BIGINT,
        INTEGER,
        LONGTEXT,
        MEDIUMTEXT,
        SET,
        SMALLINT,
        TINYINT,
        TINYTEXT,
        VARCHAR,
        )
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata

class PredictedInteractor(Base):
    __tablename__ = 'predicted_interactor'

    predicted_interactor_id = Column(INTEGER(11), primary_key=True, nullable=False, autoincrement=True)
    curated_interactor_id = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), nullable=False, index=True)
    interactor_type = Column(Enum('protein', 'gene', 'mRNA', 'synthetic'), nullable=False)
    prediction_method_id = Column(INTEGER(11), ForeignKey("prediction_method.prediction_method_id"), nullable=False, index=True)
    curies = Column(String(255))
    name = Column(String(255), nullable=False)
    molecular_structure = Column(String(10000))
    predicted_timestamp = Column(TIMESTAMP, nullable=False)
    ensembl_gene_id = Column(INTEGER(11), ForeignKey("ensembl_gene.ensembl_gene_id"), nullable=False, index=True)

    curated_interactors = relationship("CuratedInteractor", back_populates="predicted_interactors")
    prediction_methods = relationship("PredictionMethod", back_populates="predicted_interactors")
    ensembl_genes = relationship("EnsemblGene", back_populates="predicted_interactors")

    def __repr__(self):
        try:
            pi_id = self.predicted_interactor_id
            return "<PredictedInteractor(predicted_interactor_id='%d', curated_interactor_id='%d', interactor_type='%s', prediction_method_id='%d', curies='%s', name='%s', molecular_structure='%s', predicted_timestamp='%s', ensembl_gene_id='%d')>" % (
                pi_id, self.curated_interactor_id, self.interactor_type, self.prediction_method_id, self.curies, self.name, self.molecular_structure, str(self.predicted_timestamp),self.ensembl_gene_id)
        except NameError:
            return "<PredictedInteractor(predicted_interactor_id=Null-until-stored, curated_interactor_id='%d', interactor_type='%s', prediction_method_id='%d', curies='%s', name='%s', molecular_structure='%s', predicted_timestamp='%s', ensembl_gene_id='%d')>" % (
                self.curated_interactor_id, self.interactor_type, self.prediction_method_id, self.curies, self.name, self.molecular_structure, str(self.predicted_timestamp),self.ensembl_gene_id)

class CuratedInteractor(Base):
    __tablename__ = 'curated_interactor'

    curated_interactor_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    interactor_type = Column(Enum('protein', 'gene', 'mRNA', 'synthetic'), nullable=False)
    curies = Column(String(255), unique=True)   
    name = Column(String(255), nullable=True)
    molecular_structure = Column(String(10000), nullable=False)
    import_timestamp = Column(TIMESTAMP, nullable=False)
    ensembl_gene_id = Column(INTEGER(11), ForeignKey("ensembl_gene.ensembl_gene_id"), nullable=False, index=True)

    predicted_interactors = relationship("PredictedInteractor", back_populates="curated_interactors")
    interactors_1 = relationship("Interaction", primaryjoin="CuratedInteractor.curated_interactor_id == Interaction.interactor_1")
    interactors_2 = relationship("Interaction", primaryjoin="CuratedInteractor.curated_interactor_id == Interaction.interactor_2")
    ensembl_genes = relationship("EnsemblGene", back_populates="curated_interactors")

    def __repr__(self):
        try:
            ci_id = self.curated_interactor_id
            return "<CuratedInteractor(curated_interactor_id='%d', interactor_type='%s', curies='%s', name='%s', molecular_structure='%s', import_timestamp='%s', ensembl_gene_id='%d')>" % (
                ci_id, self.interactor_type, self.curies, self.name, self.molecular_structure, str(self.import_timestamp),self.ensembl_gene_id)
        except NameError:
            return "<CuratedInteractor(curated_interactor_id=Null-until-stored, interactor_type='%s', curies='%s', name='%s', molecular_structure='%s', import_timestamp='%s', ensembl_gene_id='%d')>" % (
                self.interactor_type, self.curies, self.name, self.molecular_structure, str(self.import_timestamp),self.ensembl_gene_id)

class EnsemblGene(Base):
    __tablename__ = 'ensembl_gene'

    ensembl_gene_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    species_id = Column(INTEGER(11), ForeignKey("species.species_id"), nullable=False, index=True)
    ensembl_stable_id = Column(String(255))
    import_time_stamp = Column(TIMESTAMP, nullable=False)

    predicted_interactors = relationship("PredictedInteractor", back_populates="ensembl_genes")
    curated_interactors = relationship("CuratedInteractor", back_populates="ensembl_genes")
    species_ids_r = relationship("Species", back_populates="ensembl_genes")

    def __repr__(self):
        try:
            eg_id = self.gene_id
            return "<EnsemblGene(gene_id='%d', ensembl_stable_id='%s', species_id='%d', import_time_stamp='%s')>" % (
                eg_id, self.ensembl_stable_id, self.species_id, str(self.import_time_stamp))
        except NameError:
            return "<EnsemblGene(gene_id=Null-until-stored, ensembl_stable_id='%s', species_id='%d', import_time_stamp='%s')>" % (
                self.ensembl_stable_id, self.species_id, str(self.import_time_stamp))

class Interaction(Base):
    __tablename__ = 'interaction'
    __table_args__ = (
        Index('interaction_int1_int2_doi_db', 'interactor_1', 'interactor_2', 'doi', 'source_db_id', unique=True),
    )

    interaction_id = Column(INTEGER(11), primary_key=True, nullable=False, autoincrement=True)
    interactor_1 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), nullable=False, index=True)
    interactor_2 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), nullable=False, index=True)
    doi = Column(String(255), nullable=False)   
    source_db_id = Column(INTEGER(11), ForeignKey("source_db.source_db_id"), nullable=False, index=True)
    import_timestamp = Column(TIMESTAMP, nullable=False)
    key_value_id = Column(INTEGER(11), ForeignKey("key_value_pair.key_value_id"), nullable=False, index=True)

    source_dbs = relationship("SourceDb", back_populates="interactions")
    key_value_pairs = relationship("KeyValuePair", back_populates="interactions")
    
    def __repr__(self):
        try:
            i_id = self.interaction_id
            return "<Interaction(interaction_id='%d', interactor_1='%s', interactor_2='%s', doi='%s', source_db_id='%d', import_timestamp='%s', meta_id='%d')>" % (
                i_id, self.interactor_1, self. interactor_2, self.doi, self.source_db_id, str(self.import_timestamp), self.meta_id)
        except NameError:
            return "<Interaction(interaction_id=Null-until-stored, interactor_1='%s', interactor_2='%s', doi='%s', source_db_id='%d', import_timestamp='%s', meta_id='%d')>" % (
                self.interaction_id, self.interactor_1, self. interactor_2, self.doi, self.source_db_id, str(self.import_timestamp), self.meta_id)

class KeyValuePair(Base):
    __tablename__ = 'key_value_pair'
    __table_args__ = (
        Index('key_value_pair_metakey_val', 'meta_key_id', 'value', unique=True),
    )

    key_value_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    meta_key_id = Column(INTEGER(11), ForeignKey("meta_key.meta_key_id"), nullable=False,  index=True)
    value = Column(String(255), nullable=False)
    ontology_term_id = Column(INTEGER(11), ForeignKey("ontology_term.ontology_term_id"), index=True)

    interactions = relationship("Interaction", back_populates="key_value_pairs")
    meta_keys = relationship("MetaKey", back_populates="key_value_pairs")
    ontology_terms = relationship("OntologyTerm", back_populates="key_value_pairs")

    def __repr__(self):
        try:    
            kv_id = self.key_value_id
            return "<KeyValuePair(key_value_id='%d', meta_key_id='%d', value='%s', ontology_term_id='%d')>" % (
                kv_id, self.meta_key_id, self.value, self.ontology_term_id)
        except TypeError:
            return "<KeyValuePair(key_value_id=Null-until-stored, meta_key_id='%d', value='%s', ontology_term_id='%d')>" % (
                self.meta_key_id, self.value, self.ontology_term_id)

class MetaKey(Base):
    __tablename__ = 'meta_key'

    meta_key_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False, unique=True)
    description = Column(String(255), nullable=False, unique=True)

    key_value_pairs = relationship("KeyValuePair", back_populates="meta_keys")

    def __repr__(self):
        try:    
            k_id = self.meta_key_id
            return "<MetaKey(meta_key_id='%d', name='%s', description='%s')>" % (
                k_id, self.key_name, self.key_description)
        except NameError:
            return "<MetaKey(meta_key_id=Null-until-stored, name='%s', description='%s')>" % (
                self.key_name, self.key_description)


class OntologyTerm(Base):
    __tablename__ = 'ontology_term'

    ontology_term_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    ontology_id = Column(INTEGER(11), ForeignKey("ontology.ontology_id"), nullable=False, index=True)
    accession = Column(String(255), nullable=False, unique=True)
    description = Column(String(255), nullable=False, unique=True)

    ontologies = relationship("Ontology", back_populates="ontology_terms")
    key_value_pairs = relationship("KeyValuePair", back_populates="ontology_terms")

    def __repr__(self):
        try:    
            ot_id = self.ontology_term_id
            return "<Key(ontology_term_id='%d', ontology_id='%d', accession='%s', description='%s')>" % (
                ot_id, self.ontology_id, self.accession, self.description)
        except NameError:
            return "<Key(ontology_term_id=Null-until-stored, ontology_id='%d', accession='%s', description='%s')>" % (
                self.ontology_id, self.accession, self.description)


class Ontology(Base):
    __tablename__ = 'ontology'

    ontology_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False, unique=True)
    description = Column(String(255), nullable=False, unique=True)

    ontology_terms = relationship("OntologyTerm", back_populates="ontologies")

    def __repr__(self):
        try:
            o_id = self.ontology_id
            return "<Ontology(ontology_id='%d', ontology_name='%s', ontology_description='%s')>" % (
                o_id, self.ontology_name, self.ontology_description)
        except NameError:
            return "<Ontology(ontology_id=Null-until-stored, ontology_name='%s', ontology_description='%s')>" % (
                self.ontology_name, self.ontology_description)

class PredictionMethod(Base):
    __tablename__ = 'prediction_method'

    prediction_method_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False)
    parameters = Column(String(255), nullable=False)

    predicted_interactors = relationship("PredictedInteractor", back_populates="prediction_methods")

    def __repr__(self):
        try:
            pm_id = self.prediction_method_id
            return "<PredictionMethod(prediction_method_id='%d', prediction_method_name='%s', prediction_method_values='%s')>" % (
                pm_id, self.prediction_method_name, self.prediction_method_values)
        except NameError:
            return "<PredictionMethod(prediction_method_id=Null-until-stored, prediction_method_name='%s', prediction_method_values='%s')>" % (
                self.prediction_method_name, self.prediction_method_values)
            

class SourceDb(Base):
    __tablename__ = 'source_db'

    source_db_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    label = Column(String(255), nullable=False)
    external_db = Column(String(255), nullable=False)

    interactions= relationship("Interaction", back_populates="source_dbs")

    def __repr__(self):
        try:
            sdb_id = self.source_db_id
            return "<SourceDb(source_db_id='%d', label='%s', external_db='%s')>" % (
                sdb_id, self.label, self.external_db)
        except NameError:
            return "<SourceDb(source_db_id=Null-until-stored, label='%s', external_db='%s')>" % (
                self.label, self.external_db)

class Species(Base):
    __tablename__ = 'species'

    species_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    ensembl_division = Column(String(255), nullable=False)
    production_name = Column(String(255), nullable=False)
    taxon_id = Column(INTEGER(11), nullable=False, unique=True)

    ensembl_genes = relationship("EnsemblGene", back_populates="species_ids_r")

    def __repr__(self):
        try:
            s_id = self.species_id
            return "<Species(species_id='%d', ensembl_division='%s', production_name='%s', taxon_id='%d')>" % (
                s_id, self.ensembl_division, self.production_name, self.taxon_id)
        except NameError:
            return "<Species(species_id=Null-until-stored, ensembl_division='%s', production_name='%s', taxon_id='%d')>" % (
                self.ensembl_division, self.production_name, self.taxon_id)
