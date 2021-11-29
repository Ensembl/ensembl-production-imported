DROP DATABASE P2PI_24_Nov;
CREATE DATABASE IF NOT EXISTS P2PI_24_Nov;

use P2PI_24_Nov;

CREATE TABLE `interaction` (
	`interaction_id` INT NOT NULL AUTO_INCREMENT,
	`interactor_1` INT NOT NULL,
	`interactor_2` INT NOT NULL,
	`doi` VARCHAR(255) NOT NULL,
	`source_db_id` INT NOT NULL,
	`import_timestamp` TIMESTAMP NOT NULL,
	`key_value_id` INT NOT NULL,
	PRIMARY KEY (`interaction_id`),
	UNIQUE KEY `interaction_int1_int2_doi_db` (`interactor_1`,`interactor_2`,`doi`, `source_db_id`)
);

CREATE TABLE `predicted_interactor` (
	`predicted_interaction_id` INT NOT NULL AUTO_INCREMENT,
	`curated_interactor_id` INT NOT NULL,
	`interactor_type` ENUM('protein','gene','mRNA','synthetic') NOT NULL,
	`prediction_method_id` INT NOT NULL,
	`curies` VARCHAR(255),
	`name` VARCHAR(255) NOT NULL,
	`molecular_structure` VARCHAR(10000) ,
	`predicted_timestamp` TIMESTAMP NOT NULL,
	`ensembl_gene_id` INT NOT NULL,
	PRIMARY KEY (`predicted_interaction_id`)
);

CREATE TABLE `ensembl_gene` (
	`ensembl_gene_id` INT NOT NULL AUTO_INCREMENT,
	`species_id` INT NOT NULL,
	`ensembl_stable_id` VARCHAR(255),
	`import_time_stamp` TIMESTAMP NOT NULL,
	PRIMARY KEY (`ensembl_gene_id`)
);

CREATE TABLE `curated_interactor` (
	`curated_interactor_id` INT NOT NULL AUTO_INCREMENT,
	`interactor_type` ENUM('protein','gene','mRNA','synthetic') NOT NULL,
	`curies` VARCHAR(255) UNIQUE,
	`name` VARCHAR(255),
	`molecular_structure` VARCHAR(10000) ,
	`import_timestamp` TIMESTAMP NOT NULL,
	`ensembl_gene_id` INT NOT NULL,
	PRIMARY KEY (`curated_interactor_id`)
);

CREATE TABLE `species` (
	`species_id` INT NOT NULL AUTO_INCREMENT,
	`ensembl_division` varchar(255) NOT NULL,
	`production_name` varchar(255) NOT NULL,
	`taxon_id` INT NOT NULL UNIQUE,
	PRIMARY KEY (`species_id`)
);

CREATE TABLE `key_value_pair` (
	`key_value_id` INT NOT NULL AUTO_INCREMENT,
	`key_id` INT NOT NULL,
	`value` varchar(255) NOT NULL,
	`ontology_term_id` INT,
	PRIMARY KEY (`key_value_id`)
);

CREATE TABLE `ontology` (
	`ontology_id` INT NOT NULL AUTO_INCREMENT,
	`name` VARCHAR(255) NOT NULL UNIQUE,
	`description` varchar(255) NOT NULL UNIQUE,
	PRIMARY KEY (`ontology_id`)
);

CREATE TABLE `source_db` (
	`source_db_id` INT NOT NULL AUTO_INCREMENT,
	`label` varchar(255) NOT NULL,
	`external_db` VARCHAR(255) NOT NULL,
	PRIMARY KEY (`source_db_id`)
);

CREATE TABLE `key` (
	`key_id` INT NOT NULL AUTO_INCREMENT,
	`name` VARCHAR(255) NOT NULL UNIQUE,
	`description` VARCHAR(255) NOT NULL UNIQUE,
	PRIMARY KEY (`key_id`)
);

CREATE TABLE `prediction_method` (
	`prediction_method_id` INT NOT NULL AUTO_INCREMENT,
	`name` VARCHAR(255) NOT NULL,
	`parameters` VARCHAR(255) NOT NULL,
	PRIMARY KEY (`prediction_method_id`)
);

CREATE TABLE `ontology_term` (
	`ontology_term_id` INT NOT NULL AUTO_INCREMENT,
	`ontology_id` INT NOT NULL,
	`accession` VARCHAR(255) NOT NULL UNIQUE,
	`description` VARCHAR(255) NOT NULL UNIQUE,
	PRIMARY KEY (`ontology_term_id`)
);

ALTER TABLE `interaction` ADD CONSTRAINT `interaction_fk0` FOREIGN KEY (`interactor_1`) REFERENCES `curated_interactor`(`curated_interactor_id`);

ALTER TABLE `interaction` ADD CONSTRAINT `interaction_fk1` FOREIGN KEY (`interactor_2`) REFERENCES `curated_interactor`(`curated_interactor_id`);

ALTER TABLE `interaction` ADD CONSTRAINT `interaction_fk2` FOREIGN KEY (`source_db_id`) REFERENCES `source_db`(`source_db_id`);

ALTER TABLE `interaction` ADD CONSTRAINT `interaction_fk3` FOREIGN KEY (`key_value_id`) REFERENCES `key_value_pair`(`key_value_id`);

ALTER TABLE `predicted_interactor` ADD CONSTRAINT `predicted_interactor_fk0` FOREIGN KEY (`curated_interactor_id`) REFERENCES `curated_interactor`(`curated_interactor_id`);

ALTER TABLE `predicted_interactor` ADD CONSTRAINT `predicted_interactor_fk1` FOREIGN KEY (`prediction_method_id`) REFERENCES `prediction_method`(`prediction_method_id`);

ALTER TABLE `predicted_interactor` ADD CONSTRAINT `predicted_interactor_fk2` FOREIGN KEY (`ensembl_gene_id`) REFERENCES `ensembl_gene`(`ensembl_gene_id`);

ALTER TABLE `ensembl_gene` ADD CONSTRAINT `ensembl_gene_fk0` FOREIGN KEY (`species_id`) REFERENCES `species`(`species_id`);

ALTER TABLE `curated_interactor` ADD CONSTRAINT `curated_interactor_fk0` FOREIGN KEY (`ensembl_gene_id`) REFERENCES `ensembl_gene`(`ensembl_gene_id`);

ALTER TABLE `key_value_pair` ADD CONSTRAINT `key_value_pair_fk0` FOREIGN KEY (`key_id`) REFERENCES `key`(`key_id`);

ALTER TABLE `key_value_pair` ADD CONSTRAINT `key_value_pair_fk1` FOREIGN KEY (`ontology_term_id`) REFERENCES `ontology_term`(`ontology_term_id`);

ALTER TABLE `ontology_term` ADD CONSTRAINT `ontology_term_fk0` FOREIGN KEY (`ontology_id`) REFERENCES `ontology`(`ontology_id`);




