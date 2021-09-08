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
                self.gene_id, self.ensembl_stable_id, self.species_id, str(self.import_timestamp))

class Interaction(Base):
    __tablename__ = 'interaction'

    interaction_id = Column(INTEGER(11), nullable=False)
    interactor_1 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), primary_key=True, nullable=False, index=True)
    interactor_2 = Column(INTEGER(11), ForeignKey("curated_interactor.curated_interactor_id"), primary_key=True, nullable=False, index=True)
    doi = Column(String(255), primary_key=True, nullable=False)
    source_db_id = Column(INTEGER(11), ForeignKey("source_db.source_db_id"), primary_key=True, nullable=False, index=True)
    import_timestamp = Column(TIMESTAMP, nullable=False)
    interaction_meta_id = Column(INTEGER(11), ForeignKey("interaction_meta.interaction_meta_id"), nullable=False, index=True)

    #curated_interactors_1 = relationship("CuratedInteractor", back_populates="interactors_1")
    #curated_interactors_2 = relationship("CuratedInteractor", back_populates="interactors_2")
    source_dbs = relationship("SourceDb", back_populates="interactions")
    interaction_metas = relationship("InteractionMeta", back_populates="interactions")

    def __repr__(self):
        return "<Interaction(interaction_id='%d', interactor_1='%s', interactor_2='%s', doi='%s', source_db_id='%d', import_timestamp='%s', interaction_meta_id='%d')>" % (
                self.interaction_id, self.interactor_1, self. interactor_2, self.doi, self.source_db_id, str(self.import_timestamp), self.interaction_meta_id)

class InteractionMeta(Base):
    __tablename__ = 'interaction_meta'

    interaction_meta_id = Column(INTEGER(11), primary_key=True, nullable=False)
    interaction_meta_type_id = Column(String(255), ForeignKey("interaction_meta_type.interaction_meta_type_id"), primary_key=True, nullable=False, index=True)
    reference_ontology_id = Column(INTEGER(11), ForeignKey("reference_ontology.reference_ontology_id"), nullable=False, index=True)
    ontology_accession = Column(String(255), nullable=False)
    value = Column(Float)
    ontology_decription = Column(String(255), nullable=False)

    interactions = relationship("Interaction", back_populates="interaction_metas")
    interaction_meta_types = relationship("InteractionMetaType", back_populates="interaction_metas")
    reference_ontologies = relationship("ReferenceOntology", back_populates="interaction_metas")

    def __repr__(self):
                return "<InteractionMeta(interaction_meta_id='%d', interaction_meta_type_id='%s', reference_ontology_id='%d', ontology_accession='%s', value='%d', ontology_decription='%s')>" % (
                        self.interaction_meta_id, self.interaction_meta_type_id, self.reference_ontology_id, self.ontology_accession, self.value, self.ontology_decription)

class InteractionMetaType(Base):
    __tablename__ = 'interaction_meta_type'

    interaction_meta_type_id = Column(String(255), primary_key=True)
    interaction_type = Column(String(255), nullable=False, unique=True)
    interaction_type_description = Column(String(255))

    interaction_metas = relationship("InteractionMeta", back_populates="interaction_meta_types")

    def __repr__(self):
        return "<InteractionMetaType(interaction_meta_type_id='%s', interaction_type='%s', interaction_type_description='%s')>" % (
                self.interaction_meta_type_id, self.interaction_type, self.interaction_type_description)

class PredictionMethod(Base):
    __tablename__ = 'prediction_method'

    prediction_method_id = Column(INTEGER(11), primary_key=True)
    prediction_method_name = Column(String(255), nullable=False)
    prediction_method_values = Column(String(255), nullable=False, unique=True)

    candidate_interactors = relationship("CandidateInteractor", back_populates="prediction_methods")

    def __repr__(self): 
        return "<PredictionMethod(prediction_method_id='%d', prediction_method_name='%s', prediction_method_values='%s')>" % (
                self.prediction_method_id, self.prediction_method_name, self.prediction_method_values)

class ReferenceOntology(Base):
    __tablename__ = 'reference_ontology'

    reference_ontology_id = Column(INTEGER(11), primary_key=True)
    ontology_name = Column(String(255), nullable=False, unique=True)
    ontology_description = Column(String(255), nullable=False)

    interaction_metas = relationship("InteractionMeta", back_populates="reference_ontologies")

    def __repr__(self): 
        return "<ReferenceOntology(reference_ontology_id='%d', ontology_name='%s', ontology_description='%s')>" % (
            self.reference_ontology_id, self.ontology_name, self.ontology_description)

class SourceDb(Base):
    __tablename__ = 'source_db'

    source_db_id = Column(INTEGER(11), primary_key=True)
    label = Column(String(255), nullable=False)
    external_db = Column(String(255), nullable=False)

    curated_interactors = relationship("CuratedInteractor", back_populates="source_dbs")
    interactions= relationship("Interaction", back_populates="source_dbs")

    def __repr__(self): 
        return "<SourceDb(source_db_id='%d', label='%s', external_db='%s')>" % (
                self.source_db_id, self.label, self.external_db)

class Species(Base):
    __tablename__ = 'species'

    species_id = Column(INTEGER(11), primary_key=True)
    ensembl_division = Column(String(255), nullable=False)
    species_production_name = Column(String(255), nullable=False)
    species_taxon_id = Column(INTEGER(11), nullable=False, unique=True)

    ensembl_genes = relationship("EnsemblGene", back_populates="species_ids_r")

    def __repr__(self): 
        return "<Species(species_id='%d', ensembl_division='%s', species_production_name='%s', species_taxon_id='%d')>" % (
                self.species_id, self.ensembl_division, self.species_production_name, self. species_taxon_id)
