file=open(r'D:\test\genewisetest\filtered.gff','r')
record_box={}
for record in file:
    if record.split('\t')[0] not in list(record_box.keys()):
        record_box[record.split('\t')[0]] = [record.split('\t')[3], record.split('\t')[4],
                                             record.split('\t')[5], record.split('\t')[6]]
    elif record.split('\t')[0] in list(record_box.keys()):
        if record.split('\t')[3]<=record_box[record.split('\t')[0]][0]<=record.split('\t')[4] or \
                record.split('\t')[3]<=record_box[record.split('\t')[0]][1]<=record.split('\t')[4]:
            if record.split('\t')[5] >= record_box[record.split('\t')[0]][3]:
                record_box[record.split('\t')[0]] = [record.split('\t')[3], record.split('\t')[4],
                                                     record.split('\t')[5], record.split('\t')[6]]






