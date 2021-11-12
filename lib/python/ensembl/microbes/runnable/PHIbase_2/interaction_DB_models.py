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


class CandidateInteractor(Base):
    __tablename__ = 'candidate_interactor'

    candidate_interaction_id = Column(INTEGER(11), primary_key=True, nullable=False)
    curated_interactor_id = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), nullable=False, index=True)
    interactor_type = Column(String(255), nullable=False)
    prediction_method_id = Column(INTEGER(11),ForeignKey("prediction_method.prediction_method_id"), nullable=False, index=True)
    curies = Column(String(255))
    name = Column(String(255), nullable=False)
    molecular_structure = Column(String(10000), nullable=False)
    predicted_timestamp = Column(TIMESTAMP, nullable=False)
    ensembl_gene_id = Column(INTEGER(11),ForeignKey("ensembl_gene.gene_id"), nullable=False, index=True)

    curated_interactors = relationship("CuratedInteractor", back_populates="candidate_interactors")
    prediction_methods = relationship("PredictionMethod", back_populates="candidate_interactors")
    ensembl_genes = relationship("EnsemblGene", back_populates="candidate_interactors")

    def __repr__(self):
        return "<CandidateInteractor(candidate_interaction_id='%d', curated_interactor_id='%d', interactor_type='%s', prediction_method_id='%d', curies='%s', name='%s', molecular_structure='%s', predicted_timestamp='%s', ensembl_gene_id='%d')>" % (
                self.candidate_interaction_id, self.curated_interactor_id, self.interactor_type, self.prediction_method_id, self.curies, self.name, self.molecular_structure, str(self.predicted_timestamp),self.ensembl_gene_id)

class CuratedInteractor(Base):
    __tablename__ = 'curated_interactor'

    curated_interactor_id = Column(INTEGER(11), primary_key=True)
    interactor_type = Column(String(255), nullable=False)
    curies = Column(String(255))
    name = Column(String(255), nullable=False)
    molecular_structure = Column(String(10000), nullable=False)
    import_timestamp = Column(TIMESTAMP, nullable=False)
    source_db_id = Column(INTEGER(11), ForeignKey("source_db.source_db_id"), nullable=False, index=True)
    ensembl_gene_id = Column(INTEGER(11), ForeignKey("ensembl_gene.gene_id"), nullable=False, index=True)

    candidate_interactors = relationship("CandidateInteractor", back_populates="curated_interactors")
    interactors_1 = relationship("Interaction", primaryjoin="CuratedInteractor.curated_interactor_id == Interaction.interactor_1")
    interactors_2 = relationship("Interaction", primaryjoin="CuratedInteractor.curated_interactor_id == Interaction.interactor_2")
    ensembl_genes = relationship("EnsemblGene", back_populates="curated_interactors")
    source_dbs = relationship("SourceDb", back_populates="curated_interactors")

    def __repr__(self):
        return "<CuratedInteractor(curated_interactor_id='%d', interactor_type='%s', curies='%s', name='%s', molecular_structure='%s', import_timestamp='%s', source_db_id='%d', ensembl_gene_id='%d')>" % (
                self.curated_interactor_id, self.interactor_type, self.curies, self.name, self.molecular_structure, str(self.import_timestamp),self.source_db_id, self.ensembl_gene_id)

class EnsemblGene(Base):
    __tablename__ = 'ensembl_gene'

    gene_id = Column(INTEGER(11), primary_key=True)
    ensembl_stable_id = Column(String(255))
    species_id = Column(INTEGER(11), ForeignKey("species.species_id"), nullable=False, index=True)
    import_time_stamp = Column(TIMESTAMP, nullable=False)

    candidate_interactors = relationship("CandidateInteractor", back_populates="ensembl_genes")
    curated_interactors = relationship("CuratedInteractor", back_populates="ensembl_genes")
    species_ids_r = relationship("Species", back_populates="ensembl_genes")

    def __repr__(self):
        return "<EnsemblGene(gene_id='%d', ensembl_stable_id='%s', species_id='%d', import_time_stamp='%s')>" % (
                self.gene_id, self.ensembl_stable_id, self.species_id, str(self.import_time_stamp))

class Interaction(Base):
    __tablename__ = 'interaction'

    interaction_id = Column(INTEGER(11), nullable=False)
    interactor_1 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), primary_key=True, nullable=False, index=True)
    interactor_2 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), primary_key=True, nullable=False, index=True)
    doi = Column(String(255), primary_key=True, nullable=False)
    source_db_id = Column(INTEGER(11), ForeignKey("source_db.source_db_id"), primary_key=True, nullable=False, index=True)
    import_timestamp = Column(TIMESTAMP, nullable=False)
    interaction_meta_id = Column(INTEGER(11), ForeignKey("meta_value.meta_value_id"), nullable=False, index=True)

    #curated_interactors_1 = relationship("CuratedInteractor", back_populates="interactors_1")
    #curated_interactors_2 = relationship("CuratedInteractor", back_populates="interactors_2")
    source_dbs = relationship("SourceDb", back_populates="interactions")
    meta_values = relationship("MetaValue", back_populates="interactions")

    def __repr__(self):
        return "<Interaction(interaction_id='%d', interactor_1='%s', interactor_2='%s', doi='%s', source_db_id='%d', import_timestamp='%s', interaction_meta_id='%d')>" % (
                self.interaction_id, self.interactor_1, self. interactor_2, self.doi, self.source_db_id, str(self.import_timestamp), self.interaction_meta_id)

class MetaValue(Base):
    __tablename__ = 'meta_value'

    meta_value_id = Column(INTEGER(11), primary_key=True, nullable=False)
    meta_key_id = Column(String(255), ForeignKey("meta_key.meta_key_id"),  primary_key=True, nullable=False, index=True)
    value = Column(String(255), nullable=False)
    ontology_accession = Column(String(255))
    reference_ontology_id = Column(INTEGER(11), ForeignKey("meta_ontology.meta_ontology_id"), index=True)
    float_value = Column(Float)

    interactions = relationship("Interaction", back_populates="meta_values")
    meta_keys = relationship("MetaKey", back_populates="meta_values")
    meta_ontologies = relationship("MetaOntology", back_populates="meta_values")

    def __repr__(self):
                return "<MetaValue(meta_value_id='%d', meta_key_id='%s', value='%s', reference_ontology_id='%d', ontology_accession='%s', float_value='%d')>" % (
                        self.meta_value_id, self.meta_key_id, self.value, self.reference_ontology_id, self.ontology_accession, self.float_value)

class MetaKey(Base):
    __tablename__ = 'meta_key'

    meta_key_id = Column(String(255), primary_key=True)
    key_name = Column(String(255), nullable=False, unique=True)
    key_description = Column(String(255))

    meta_values = relationship("MetaValue", back_populates="meta_keys")

    def __repr__(self):
        return "<MetaKey(meta_key_id='%s', key_name='%s', key_description='%s')>" % (
                self.meta_key_id, self.key_name, self.key_description)

class PredictionMethod(Base):
    __tablename__ = 'prediction_method'

    prediction_method_id = Column(INTEGER(11), primary_key=True)
    prediction_method_name = Column(String(255), nullable=False)
    prediction_method_values = Column(String(255), nullable=False, unique=True)

    candidate_interactors = relationship("CandidateInteractor", back_populates="prediction_methods")

    def __repr__(self): 
        return "<PredictionMethod(prediction_method_id='%d', prediction_method_name='%s', prediction_method_values='%s')>" % (
                self.prediction_method_id, self.prediction_method_name, self.prediction_method_values)

class MetaOntology(Base):
    __tablename__ = 'meta_ontology'

    meta_ontology_id = Column(INTEGER(11), primary_key=True)
    ontology_name = Column(String(255), nullable=False, unique=True)
    ontology_description = Column(String(255), nullable=False)

    meta_values = relationship("MetaValue", back_populates="meta_ontologies")

    def __repr__(self): 
        return "<MetaOntology(meta_ontology_id='%d', ontology_name='%s', ontology_description='%s')>" % (
            self.reference_ontology_id, self.ontology_name, self.ontology_description)

class SourceDb(Base):
    __tablename__ = 'source_db'

    source_db_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    label = Column(String(255), nullable=False)
    external_db = Column(String(255), nullable=False)

    curated_interactors = relationship("CuratedInteractor", back_populates="source_dbs")
    interactions= relationship("Interaction", back_populates="source_dbs")

    def __repr__(self): 
        return "<SourceDb(label='%s', external_db='%s')>" % (
                self.label, self.external_db)

class Species(Base):
    __tablename__ = 'species'

    species_id = Column(INTEGER(11), primary_key=True, autoincrement=True)
    ensembl_division = Column(String(255), nullable=False)
    species_production_name = Column(String(255), nullable=False)
    species_taxon_id = Column(INTEGER(11), nullable=False, unique=True)

    ensembl_genes = relationship("EnsemblGene", back_populates="species_ids_r")

    def __repr__(self): 
        return "<Species(ensembl_division='%s', species_production_name='%s', species_taxon_id='%d')>" % (
                self.ensembl_division, self.species_production_name, self. species_taxon_id)
