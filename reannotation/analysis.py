from collections import Counter


def interpro_accessions_frequently_missed_by_all_tools(braker3_results, helixer_results, anno_results):
    for acc in set(anno_results["l3_more_frequent"].keys()).intersection(braker3_results["l3_more_frequent"].keys()).intersection(helixer_results["l3_more_frequent"].keys()):
        acc_product = braker3_results["acc_product"]
        print(f"{acc}: {acc_product[acc]}")
        acc_list = braker3_results["acc_list3"]
        freq = Counter(acc_list)[acc]
        l3_more_frequent = braker3_results["l3_more_frequent"]
        print(f"{round(l3_more_frequent[acc], 2)} times more likely than with BRAKER3 ({freq} occurrences, {round(freq/l3_more_frequent[acc])} expected)")

        acc_list = helixer_results["acc_list3"]
        freq = Counter(acc_list)[acc]
        l3_more_frequent = helixer_results["l3_more_frequent"]
        print(f"{round(l3_more_frequent[acc], 2)} times more likely than with Helixer ({freq} occurrences, {round(freq/l3_more_frequent[acc])} expected)")

        acc_list = anno_results["acc_list3"]
        freq = Counter(acc_list)[acc]
        l3_more_frequent = anno_results["l3_more_frequent"]
        print(f"{round(l3_more_frequent[acc], 2)} times more likely than with Anno ({freq} occurrences, {round(freq/l3_more_frequent[acc])} expected)")
        print("~~~~~~~~~~")
