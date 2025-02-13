from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os


driver = webdriver.Edge()

url = 'https://zhanggroup.org/ATPbind/'
driver.get(url)

email = ''  # 邮箱 (email)
folder_path = r''  # 文件夹路径 (folder path)

for file_name in os.listdir(folder_path):
    file_path = os.path.join(folder_path, file_name)

    if not os.path.isfile(file_path):
        continue

    try:
        # 上传按钮 (upload button)
        file_input = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'str_file'))
        )
        file_input.send_keys(file_path)

        # 文件名 (file name)
        file_name_input = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'TARGET-NAME'))
        )
        file_name_input.clear()
        file_name_input.send_keys(file_name)

        # 邮箱 (email)
        email_input = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'REPLY-E-MAIL'))
        )
        email_input.clear()
        email_input.send_keys(email)

        # 提交 (submit)
        submit_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@type="submit" and @value="Run ATPbind"]'))
        )
        submit_button.click()

        time.sleep(5)


        driver.get(url)

    except Exception as e:
        print(f"提交文件 {file_name} 时出错: {e}")  # 提交文件 (submit file) 时出错 (error occurred)
        continue


driver.quit()
