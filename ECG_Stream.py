import streamlit as st
import numpy as np
import cv2
from tensorflow import keras
from tensorflow.keras.applications import EfficientNetB0,EfficientNetB3

# Load the model
model = keras.models.load_model("ECG_ECN.h5")

# Define label mapping
label_mapping = {
    0: 'ECG Images of Myocardial Infarction Patients (240x12=2880)',
    1: 'ECG Images of Patient that have History of MI (172x12=2064)',
    2: 'ECG Images of Patient that have abnormal heartbeat (233x12=2796)',
    3: 'Normal Person ECG Images (284x12=3408)'
}

def preprocess_image(image_bytes):
    # Convert BytesIO to numpy array
    image_array = np.asarray(bytearray(image_bytes.read()), dtype=np.uint8)

    # Decode image
    image = cv2.imdecode(image_array, cv2.IMREAD_COLOR)
    
    # Preprocess the image
    image = image[300:1480, 150:2125]
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    image = cv2.resize(image, (224, 224))
    return image

def predict(image_bytes):
    # Preprocess the image
    image = preprocess_image(image_bytes)
    image = np.expand_dims(image, axis=0)  # Add batch dimension

    # Make prediction
    prediction = model.predict(image)
    predicted_class = np.argmax(prediction)

    return label_mapping[predicted_class], prediction[0]

def main():
    st.title("ECG Image Classifier")

    uploaded_file = st.file_uploader("Choose an ECG image...", type=["jpg", "jpeg", "png"])

    if uploaded_file is not None:
        st.image(uploaded_file, caption="Uploaded Image.", use_column_width=True)
        st.write("")
        st.write("Classifying...")

        # Make prediction
        class_prediction, class_probabilities = predict(uploaded_file)

        st.write(f"Prediction: {class_prediction}")
        st.write("Class Probabilities:")
        for i, prob in enumerate(class_probabilities):
            st.write(f"{label_mapping[i]}: {prob:.2%}")

if __name__ == "__main__":
    main()