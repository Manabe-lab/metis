import time
import sys
import os
import shutil
import glob
import re

def clear_old_directories(path, date=1):
    now = time.time()
    old_dt = now - date * 86400
    for curdir, dirs, files in os.walk(path):
        for dir in dirs:
            target = os.path.join(curdir, dir)
            if os.path.getmtime(target) < old_dt:
                print(target)
                shutil.rmtree(target, ignore_errors=True)


def clear_old_files(path, date=1):
    now = time.time()
    old_dt = now - date * 86400
    files = glob.glob(path + "/*")
    for file in files:
        epSec = os.path.getctime(file)
        if epSec < old_dt:
            os.remove(file)


# temp / res directoryを作る temp_pathを指定するとtemp_pathは作らない
# 戻り値の最後に/はつかない。
def mk_temp_dir(res_path, temp_path=None):
    if temp_path is None:
        temp_dir = "temp/" + str(round(time.time()))
        if not os.path.exists('temp'):
            os.mkdir('temp')
        else:
            clear_old_directories("temp")
            clear_old_files("temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    else:
        temp_dir = temp_path
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    res_dir = temp_dir + '/' +  res_path
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
    return temp_dir, res_dir




def check_species(lst):
    if not lst:
        return False
    list_str = ''.join(lst[:50])
    uppercase_chars = sum(1 for char in ''.join(lst) if char.isupper())
    if uppercase_chars/len(list_str) > 0.6:
    	return 'human'
    else:
    	return 'mouse'


def check_species_index(lst):
    if not lst:
        return False
    list_str = ''.join(lst[:50])
    uppercase_chars = sum(1 for char in ''.join(lst) if char.isupper())
    if uppercase_chars/len(list_str) > 0.6:
    	return 1
    else:
    	return 0


def make_r_compatible_column_names(df):
    def clean_column_name(name):
        # 先頭が数字の場合、X_を追加
        if name[0].isdigit():
            name = 'X_' + name
        
        # スペースをアンダースコアに置換
        name = name.replace(' ', '_')
        
        # +をplusに、-をminusに置換
        name = name.replace('+', 'Plus')
        
        # 特殊文字を削除（アンダースコアは保持）
        name = re.sub(r'[^a-zA-Z0-9_]', '', name)
        
        # 名前が空になった場合、X_を使用
        if not name:
            name = 'X_'
        
        return name

    # すべてのカラム名に対して clean_column_name を適用
    new_columns = {col: clean_column_name(col) for col in df.columns}
    
    # 重複するカラム名を処理
    seen = set()
    for old_name, new_name in new_columns.items():
        if new_name in seen:
            i = 1
            while f"{new_name}_{i}" in seen:
                i += 1
            new_columns[old_name] = f"{new_name}_{i}"
        seen.add(new_columns[old_name])

    # データフレームのカラム名を変更
    df.rename(columns=new_columns, inplace=True)
    
    return df

def remove_after_space(i):
    # まず括弧とその中身を削除
    pattern = r'\s*[\(\[\{\<][^\)\]\}\>]*[\)\]\}\>]'
    while re.search(pattern, i):
        i = re.sub(pattern, '', i)
    # 残りのスペースを_に変換
    return re.sub(r'\s+', '_', i.strip())

def remove_sample_num(i):
    # 末尾の数字を除去
    i = re.sub('\d+$', '', i)
    # [_-]の後に続く数字を除去
    i = re.sub('[_-][^-_]*\d+$', '', i)
    # 末尾の _ または - を除去（複数回繰り返す）
    while i.endswith('_') or i.endswith('-'):
        i = i.rstrip('_-')
    return i

def remove_underscore(i):
    i = re.sub('[_-][^-_]*.+$', '', i)
    return i

def select_underscore(i):
    m = re.match(r'.+[_-]([^-_]*.+)$', i)
    if m is not None:
        return m.group(1)
    else:
        return '1'


if __name__ == '__main__':
	print("original functions")