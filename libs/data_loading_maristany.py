import json
import logging

import setup_logging

setup_logging.setup_logging()

logger = logging.getLogger(__name__)


def _load_data(file_path):
    with open(file_path) as f:
        return json.load(f)


def _glue_seq_data(seq_entries):
    seq_base = ""
    seq_muts = {}
    errors = []
    for i, seq_entry in enumerate(seq_entries):
        seq_base += seq_entry["seq_base"]
        for mut_name, mut_seq in seq_entry["seq_mut"].items():
            if mut_name not in seq_muts and i == 0:
                seq_muts[mut_name] = mut_seq
            elif mut_name not in seq_muts and i > 0:
                errors.append(f"Glue error, seq: {seq_base}, mut_name: {mut_name} was already added, block")
            else:
                seq_muts[mut_name] += mut_seq

    return seq_base, seq_muts, errors


def _validate_seqs_length(seq_base, seq_muts):
    base_length = len(seq_base)
    errors = []
    for mut_name, mut_data in seq_muts.items():
        if len(mut_data) != base_length:
            errors.append(f'Mutation lenght: {len(mut_data)} != base length: {base_length}')
    return errors


def _validate_seqs_data(seq_base, seq_muts):
    base_length = len(seq_base)
    # format - "mut_name, line_numer, error_name"
    errors = []
    for i, (mut_name, mut_data) in enumerate(seq_muts.items()):
        if len(mut_data) != base_length:
            errors[mut_name] = "wrong lenght"
        minus_started = False
        plus_started = False
        counter = ""
        value = ""

        minus_count = 0
        minus_value = None

        plus_count = 0
        plus_value = None

        for j, letter in enumerate(mut_name):
            if j == 0 and letter != "-":
                errors.append({"mut_name": mut_name,
                               "line_n": i,
                               "error_name": "mut_name_error",
                               "error_message": "unexpected minus sign"})
                continue
            elif letter == "-":
                minus_started = True
            elif letter == "+":
                minus_started = False
                plus_started = True

                if counter == "":
                    counter = "1"
                minus_count = int(counter)
                minus_value = value

                counter = ""
                value = ""

            elif minus_started or plus_started:
                if letter.isdigit():
                    counter += letter
                elif letter.isalpha():
                    value += letter

        if plus_started:
            if counter == "":
                counter = "1"
            plus_count = int(counter)
            plus_value = value
        else:
            minus_count = int(counter)
            minus_value = value

        if not minus_count or not minus_value:
            errors.append({"mut_name": mut_name,
                           "line_n": i,
                           "error_name": "mut_name_error",
                           "error_message": "no minus value or count"})
            continue

        if plus_value:
            # print(plus_count)
            # print(minus_count)
            if plus_count != minus_count:
                errors.append({"mut_name": mut_name,
                               "line_n": i,
                               "error_name": "mut_name_error",
                               "error_message": "plus_count != minus_count"})
                continue

            if mut_data.count(plus_value) != plus_count:
                errors.append({"mut_name": mut_name,
                               "line_n": i,
                               "error_name": "mut_name_error",
                               "error_message": f"wrong data - plus count ({plus_value}) - parsed:"
                                                f" {plus_count} - "
                                                f"real-value:"
                                                f" {mut_data.count(plus_value)}"})

        total_sub = 0
        for letter in mut_data:
            if letter.isalpha():
                total_sub += 1
        if total_sub != minus_count:
            if not plus_value or total_sub != mut_data.count(plus_value):
                errors.append({"mut_name": mut_name,
                               "line_n": i,
                               "error_name": "mut_name_error",
                               "error_message": f"wrong data - minus count ({minus_value}) - parsed"
                                                f"{minus_count} - "
                                                f"real-value:"
                                                f" {total_sub}"})
    return errors


def _gen_seq_from_mut(seq_base, seq_muts):
    generated_sequences = {}
    for seq_name, seq_mut in seq_muts.items():
        new_seq = ""
        for b_l, m_l in zip(seq_base, seq_mut):
            if m_l.isalpha():
                new_seq += m_l
            else:
                new_seq += b_l
        if seq_name in generated_sequences:
            raise ValueError(f"Duplicate mutation encountered - {seq_name}")
        generated_sequences[seq_name] = new_seq
    return generated_sequences


def load_and_validate_data(file_path, ignore_validation_errors=True, logging_active=True):
    data = _load_data(file_path)
    parsed_data = []
    if logging_active:
        logger.info(f"Ignoring errors - {ignore_validation_errors}")
    for gene in data.keys():
        if logging_active:
            logger.info(f'gene: {gene} - seq name: {data[gene]["seq_name"]}')
            # for entry in data[key]["seq_data"]:

        seq_base, seq_muts, glue_errors = _glue_seq_data(data[gene]["seq_data"])
        if glue_errors:
            logger.error(f'Glue error(s) encountered: Gene {gene} - {data[gene]["seq_name"]}')
            logger.error(glue_errors)
            logger.error("Breaking processing")
            raise ValueError("Glue error(s) encounter")

        length_errors = _validate_seqs_length(seq_base, seq_muts)
        if length_errors:
            logger.error(f'Length error(s) encountered: Gene {gene} - {data[gene]["seq_name"]}')
            logger.error(length_errors)
            logger.error("Breaking processing")
            raise ValueError("Length error(s) encounter")

        validation_errors = _validate_seqs_data(seq_base, seq_muts)
        if logging_active:
            logger.warning(f'Validation error(s) encountered: Gene {gene} - {data[gene]["seq_name"]}')
            logger.warning(json.dumps({"errors": validation_errors}, indent=4))
        if validation_errors and ignore_validation_errors:
            logger.info("ignoring mutations with errors")
            keys_to_del = set([error["mut_name"] for error in validation_errors])
            for key in keys_to_del:
                del seq_muts[key]
        else:
            logger.error("Breaking processing")
            raise ValueError("Validation error(s) encounter")
        seq_gen = _gen_seq_from_mut(seq_base, seq_muts)
        parsed_data.append({"gene": gene,
                            "seq_name": data[gene]["seq_name"],
                            "seq_base": seq_base,
                            "seq_gen": seq_gen})
    return parsed_data


if __name__ == "__main__":
    load_and_validate_data("../data/seqs_maristany.json")
