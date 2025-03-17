from collections import Counter

from reannotation.utils import extract_accessions_from_transcript


def interpro_accessions_frequently_missed_by_all_tools(acc_product, acc_tally, test_results, min_freq=1):
    # InterPro accessions from transcripts missed by all tools, sorted by greatest odds ratio
    print("InterPro accessions occurring with significantly higher frequency in transcripts that were missed by all tools, than in transcripts shared by at least 1 tool:")
    for acc, [stat, p] in sorted(test_results["more_frequent"].items(), key=lambda x: x[1], reverse=True):
        if acc in acc_product:
            freq = Counter(acc_tally)[acc]
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
    print()

    print("InterPro accessions that are completely missing from transcripts shared by at least 1 tool, but present in transcripts that were missed by all tools:")
    for acc, freq in Counter(acc_tally).most_common():
        if acc in test_results["not_occurring"] and acc in acc_product:
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
    print()

    print("InterPro accessions occurring as expected in transcripts that were missed by all tools.")
    for acc, [stat, p] in sorted(test_results["as_expected"].items(), key=lambda x: x[1], reverse=True):
        if acc in acc_product:
            freq = Counter(acc_tally)[acc]
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")


def interpro_accessions_frequently_missed_by_each_tool(acc_product, tools_missed_results, acc_tally_missed_tools):
    for acc in set.intersection(*[set(results["more_frequent"].keys()) for results in tools_missed_results.values()]):
        print(f"{acc}: {acc_product[acc]}")
        for tool in tools_missed_results.keys():
            acc_list = acc_tally_missed_tools[tool]
            freq = Counter(acc_list)[acc]
            more_frequent = tools_missed_results[tool]["more_frequent"]
            print(f"{round(more_frequent[acc][0], 2)} times more likely than with {tool} ({freq} occurrences, {round(freq/more_frequent[acc][0])} expected)")
        print("~~~~~~~~~~")


def interpro_accessions_in_novel_transcripts(acc_product, acc_tally, test_results, min_freq):
    # InterPro accessions from novel transcripts, sorted by greatest odds ratio
    print("InterPro accessions occurring with significantly higher frequency in novel transcripts than in shared transcripts:")
    for acc, [stat, p] in sorted(test_results["more_frequent"].items(), key=lambda x: x[1], reverse=True):
        if acc in acc_product:
            freq = Counter(acc_tally)[acc]
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
                yield acc
    print()

    print("InterPro accessions that are completely missing from shared transcripts, with high frequency in novel transcripts:")
    for acc, freq in Counter(acc_tally).most_common():
        if acc in acc_product:
            if acc in test_results["not_occurring"]:
                if freq >= min_freq:
                    print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
                    yield acc
    print()


def interpro_accessions_in_missed_transcripts(acc_product, acc_tally, acc_tally_novel, test_results, novel_test_results, min_freq):
    print("InterPro accessions that are completely missing from shared transcripts, with high frequency in missed transcripts:")
    for acc, freq in Counter(acc_tally).most_common():
        if acc in test_results["not_occurring"]:
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
    print()

    # InterPro accessions missed by automated tool, sorted by greatest odds ratio
    print("InterPro accessions occurring with significantly higher frequency in missed transcripts than in shared transcripts:")
    for acc, [stat, p] in sorted(test_results["more_frequent"].items(), key=lambda x: x[1], reverse=True):
        freq = Counter(acc_tally)[acc]
        if freq >= min_freq:
            print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
            if acc in novel_test_results["more_frequent"]:
                novel_freq = Counter(acc_tally_novel)[acc]
                novel_stat = novel_test_results["more_frequent"][acc]
                print(f"\t\t{acc} also significantly more frequent in novel transcripts ({novel_freq} occurrences, {round(novel_freq / novel_stat[0])} expected)")
    print()

    print("InterPro accessions occurring as expected in missed transcripts with high frequency:")
    for acc, [stat, p] in test_results["as_expected"].items():
        if stat < 1.5 and stat >= 0.5:
            freq = Counter(acc_tally)[acc]
            if freq >= min_freq:
                print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences)")
    print()

    print("InterPro accessions occurring less frequently in missed transcripts than expected:")
    for acc, [stat, p] in test_results["less_frequent"].items():
        freq = Counter(acc_tally)[acc]
        print(f"\t{acc}: {acc_product[acc]} ({freq} occurrences, {round(freq/stat)} expected)")
    print()


def missed_transcripts_with_significantly_more_frequent_accessions(db, missed_transcripts, acc_tally, test_results, min_freq):
    # Filter accessions with high frequency and high test statistic [print out is easy to copy/paste for Artemis].
    filt = [acc for acc, freq in Counter(acc_tally).items() if freq >= min_freq]
    filt = [acc for acc, [stat, p] in test_results["more_frequent"].items() if stat > 5 and acc in filt]

    for tran in db.all_features(featuretype="mRNA"):
        if tran.id.strip("transcript:") in missed_transcripts:
            tran_filt_accs = set(acc for acc, _ in extract_accessions_from_transcript(tran) if acc in filt)
            if tran_filt_accs:
                print(f"{tran.seqid} - {tran.id.strip('transcript:')} - {tran_filt_accs}")
